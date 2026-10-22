# -*- coding: utf-8 -*-

"""Functions for handling SMILES"""

import logging

import numpy as np

try:
    from openeye import oechem
except ImportError:
    oechem_available = False
    oechem_licensed = False
else:
    oechem_available = True
    oechem_licensed = oechem.OEChemIsLicensed()

try:
    from openeye import oeomega
except ImportError:
    oeomega_available = False
    oeomega_licensed = False
else:
    oeomega_available = True
    oeomega_licensed = oeomega.OEOmegaIsLicensed()

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdDistGeom
    from rdkit.Chem import rdmolops
except ModuleNotFoundError:
    print("Please install RDKit using conda:\n     conda install -c conda-forge rdkit")
    raise

logger = logging.getLogger(__name__)


# The random seed for RDKit's embedding of SMILES (see from_smiles)
EMBEDDING_SEED = 20261004
#: The least distance between fragments of a SMILES placed side by side, Å
FRAGMENT_GAP = 3.0

RDKit_Embedding_Error = {
    "INITIAL_COORDS": """
generation of the initial coordinates from the random distance matrix (default)
or from a set of random coordinates (when using random coordinate embedding) failed.
""",
    "FIRST_MINIMIZATION": """
the initial optimization of the atom positions using the distance-geometry force field
failed to produce a low-enough energy conformer. The check here has thresholds for both
average energy per atom and the individual atom energies.
""",
    "CHECK_TETRAHEDRAL_CENTERS": """
at least one tetrahedral C and N centers either has a volume around it which is too
small or is outside the volume defined by its neighbors
""",
    "MINIMIZE_FOURTH_DIMENSION": """
the minimization to force the values of the fourth-dimensional component of each atom
position failed
""",
    "ETK_MINIMIZATION": """
after the minimization with the ET and/or K terms, at least one atom which should have
been planar was not
""",
    "FINAL_CHIRAL_BOUNDS": """
the neighborhood of an atom with specified chirality was too distorted
(it violated distance constraints)
""",
    "FINAL_CENTER_IN_VOLUME": """
an atom with specified chirality was outside of the volume defined by its neighbors
""",
    "LINEAR_DOUBLE_BOND": """
one of the end atoms of a double bond had a linear geometry
""",
    "BAD_DOUBLE_BOND_STEREO": """
the stereochemistry of a double bond with specified stereochemistry was wrong in the
generated conformer
""",
}


def check_openeye_license():
    if not oechem_available:
        raise RuntimeError(
            "OpenEye OEChem is not installed! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )
    if not oechem_licensed:
        raise RuntimeError(
            "Cannot find a license for OpenEye OEChem! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )
    if not oeomega_available:
        raise RuntimeError(
            "OpenEye OEOmega is not installed! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )
    if not oeomega_licensed:
        raise RuntimeError(
            "Cannot find a license for OpenEye OEOmega! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )


class SMILESMixin:
    """A mixin for handling SMILES."""

    @property
    def canonical_smiles(self):
        """Return the canonical SMILES string for this object."""
        return self.to_smiles(canonical=True)

    @property
    def isomeric_smiles(self):
        """Return the canonical, isomeric SMILES string for this object."""
        return self.to_smiles(canonical=True, isomeric=True)

    @property
    def smarts(self):
        """Return the SMARTS string for this object."""
        return self.to_smarts()

    @property
    def smiles(self):
        """Return the SMILES string for this object."""
        return self.to_smiles()

    def to_smiles(
        self,
        canonical=False,
        hydrogens=False,
        isomeric=True,
        flavor="rdkit",
    ):
        """Create the SMILES string from the system.

        Parameters
        ----------
        canonical : bool = False
            Whether to create canonical SMILES
        hydrogens : bool = False
            Whether to keep H's in the SMILES string.
        isomeric : bool = True
            Whether to use isomeric SMILES
        rdkit : bool = False
            Whether to use RDKit rather than default of OpenBabel

        Returns
        -------
        str
            The SMILES string, or (SMILES, name) if the rname is requested
        """
        if flavor == "rdkit":
            mol = self.to_RDKMol()
            if isomeric:
                rdmolops.AssignStereochemistryFrom3D(mol)
            if hydrogens:
                mol2 = mol
            else:
                mol2 = AllChem.RemoveHs(mol)
            if canonical:
                smiles = Chem.MolToSmiles(mol2, isomericSmiles=isomeric)
            else:
                smiles = Chem.MolToSmiles(
                    mol2, isomericSmiles=isomeric, canonical=False
                )
        elif flavor == "openbabel":
            obConversion = openbabel.OBConversion()
            if canonical:
                obConversion.SetOutFormat("can")
            else:
                obConversion.SetOutFormat("smi")

            mol = self.to_OBMol()

            if hydrogens:
                obConversion.AddOption("h")
            if not isomeric:
                obConversion.AddOption("i")
            smiles = obConversion.WriteString(mol)
        elif flavor == "openeye":
            check_openeye_license()

            mol = self.to_OEGraphMol()
            smiles = oechem.OECreateSmiString(mol)
        else:
            raise ValueError(f"flavor of SMILES '{flavor}' not supported")

        return smiles.strip()

    def from_smiles(self, smiles, name=None, reorient=True, flavor="rdkit"):
        """Create the system from a SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string
        name : str = None
            The name of the molecule
        reorient : bool = True
            Whether to reorient to the standard orientation
        rdkit : bool = False
            Whether to use RDKit rather than default of OpenBabel

        Returns
        -------
        None
        """

        save = self.name

        if flavor == "rdkit":
            mol = Chem.rdmolfiles.MolFromSmiles(smiles)
            if mol is None:
                if Chem.rdmolfiles.MolFromSmiles(smiles, sanitize=False) is None:
                    raise ValueError(f"SMILES '{smiles}' is not valid.")
                # Valid, but RDKit refuses the valences, e.g. AsF6-
                logger.warning(
                    f"RDKit cannot handle the valences in '{smiles}', so Open Babel "
                    "builds the structure."
                )
                flavor = "openbabel"
        if flavor == "rdkit":
            mol = Chem.AddHs(mol)
            # A fixed seed: the same SMILES gives the same structure every time,
            # whatever was built before it in the process (RDKit's default
            # carries its random state from call to call), so a resumed job or a
            # parallel loop iteration builds what an uninterrupted run would.
            conformer = AllChem.EmbedMolecule(mol, randomSeed=EMBEDDING_SEED)
            if conformer == -1:
                # Check if there are small rings
                for ring in rdmolops.GetSSSR(mol):
                    if len(ring) <= 4:
                        ps = rdDistGeom.srETKDGv3()
                        break
                else:
                    ps = rdDistGeom.ETKDGv3()
                ps.trackFailures = True
                ps.randomSeed = EMBEDDING_SEED
                conformer = rdDistGeom.EmbedMolecule(mol, ps)
            if conformer == -1:
                counts = ps.GetFailureCounts()
                for i, k in enumerate(rdDistGeom.EmbedFailureCauses.names):
                    if counts[i] > 0:
                        logger.warning(
                            f" {counts[i]:5} {k} errors in embedding"
                            f"{RDKit_Embedding_Error[k]}"
                        )
                # e.g. a metal complex such as [Co(NH3)6]3+
                logger.warning(
                    f"RDKit could not embed '{smiles}', so Open Babel builds the "
                    "structure."
                )
                flavor = "openbabel"
            else:
                self.from_RDKMol(mol, properties=None)
        if flavor == "openbabel":
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("smi", "mdl")
            obConversion.AddOption("3")
            mol = openbabel.OBMol()
            obConversion.ReadString(mol, smiles)

            # Add hydrogens
            mol.AddHydrogens()

            # Get coordinates for a 3-D structure
            builder = openbabel.OBBuilder()
            builder.Build(mol)

            self.from_OBMol(mol, properties=None)
        elif flavor == "openeye":
            check_openeye_license()
            mol = oechem.OEGraphMol()
            if not oechem.OESmilesToMol(mol, smiles):
                raise ValueError(f"SMILES '{smiles}' is not valid.")
            omega = oeomega.OEOmega()
            omega.SetMaxConfs(1)
            omega.SetStrictStereo(False)
            omega.SetStrictAtomTypes(False)

            # Add hydrogens
            oechem.OEAddExplicitHydrogens(mol)

            self.from_OEMol(mol, properties=None)

        # Hypervalent and coordination centres (PF6-, SF6, metal complexes) are
        # embedded badly, so place their ligands by VSEPR; then a SMILES of several
        # fragments ('[Na+].O'), embedded on top of each other, side by side.
        self._regularize_centres()
        self._separate_fragments()

        # Rotate to standard orientation
        if reorient:
            try:
                rdkMol = self.to_RDKMol()
                rdkConf = rdkMol.GetConformers()[0]
                Chem.rdMolTransforms.CanonicalizeConformer(rdkConf, ignoreHs=False)
                self.coordinates_from_RDKMol(rdkMol)
            except Exception:
                # RDKit refuses some valences (AsF6-): the principal axes instead
                self._principal_axes_orientation()

        if name is not None:
            self.name = name
        else:
            self.name = save

    def _regularize_centres(self):
        """Place the ligands of hypervalent and coordination centres by VSEPR, and
        check that no atoms are left too close (see :mod:`molsystem.vsepr`)."""
        from . import vsepr

        n = self.n_atoms
        if n < 3:
            return
        symbols = self.atoms.symbols
        index = {atom_id: k for k, atom_id in enumerate(self.atoms.ids)}
        bonds = [
            (index[i], index[j], order)
            for i, j, order in zip(
                self.bonds.get_column_data("i"),
                self.bonds.get_column_data("j"),
                self.bonds.get_column_data("bondorder"),
            )
        ]
        charges = None
        if "formal_charge" in self.atoms:
            charges = self.atoms.get_column_data("formal_charge")
        xyz = self.atoms.get_coordinates(fractionals=False)
        xyz, changed = vsepr.regularize(symbols, xyz, bonds, charges)
        if not changed:
            return
        self.atoms.set_coordinates(xyz.tolist(), fractionals=False)
        # Fail loudly rather than return a broken structure
        molecules = self.find_molecules(as_indices=True)
        atoms = [a for m in molecules if set(m) & set(changed) for a in m]
        vsepr.check_contacts(symbols, xyz, bonds, atoms=atoms)

    def _principal_axes_orientation(self):
        """Centre the structure and align its principal axes with x, y and z."""
        xyz = np.array(self.atoms.get_coordinates(fractionals=False), dtype=float)
        xyz -= xyz.mean(axis=0)
        _, vectors = np.linalg.eigh(xyz.T @ xyz)
        xyz = xyz @ vectors[:, ::-1]
        self.atoms.set_coordinates(xyz.tolist(), fractionals=False)

    def _separate_fragments(self, gap=FRAGMENT_GAP):
        """Place the fragments of the structure side by side along x.

        RDKit and Open Babel embed the fragments of a SMILES such as '[Na+].O' or
        'O.O' each at the origin, so their atoms overlap (2026-10-21: Na 0.38 Å from
        O). Each fragment after the first is moved so that its atoms are at least
        ``gap`` Å beyond those already placed along x, centred on them in y and z: a
        start for an optimization, not an equilibrium structure.
        """
        if self.n_atoms < 2:
            return
        fragments = self.find_molecules(as_indices=True)
        if len(fragments) < 2:
            return
        xyz = np.array(self.atoms.get_coordinates(fractionals=False), dtype=float)
        placed = list(fragments[0])
        for fragment in fragments[1:]:
            fragment = list(fragment)
            shift = np.zeros(3)
            shift[0] = xyz[placed, 0].max() + gap - xyz[fragment, 0].min()
            shift[1:] = xyz[placed, 1:].mean(axis=0) - xyz[fragment, 1:].mean(axis=0)
            xyz[fragment] += shift
            placed.extend(fragment)
        self.atoms.set_coordinates(xyz.tolist(), fractionals=False)

    def to_smarts(self):
        """Generate a SMARTS string for this object.

        Returns
        -------
        str
            The SMARTS string.
        """

        generator = GenSMARTS(self)
        return generator.smarts


class GenSMARTS(object):
    """A class to generate SMARTS strings for an object.

    Parameters
    ----------
    mol_object : _Configuration, _Template, _Subset
    """

    def __init__(self, mol_object=None):
        self.smarts = ""
        self._mol_object = None
        self.symbols = []
        self.visited = []
        self.neighbors = []

        self.mol_object = mol_object

    @property
    def mol_object(self):
        """The molecular object to work with."""
        return self._mol_object

    @mol_object.setter
    def mol_object(self, value):
        self._mol_object = value

        if self._mol_object is not None:
            atoms = self.mol_object.atoms
            bonds = self.mol_object.bonds

            self.symbols = atoms.symbols

            # Find the neighbors of each atom for walking the structure
            tmp_neighbors = {i: [] for i in atoms.ids}
            for bond in bonds.bonds():
                i = bond["i"]
                j = bond["j"]
                bo = bond["bondorder"]
                tmp_neighbors[i].append((j, bo))
                tmp_neighbors[j].append((i, bo))

            # The index of each atom id in the list, 0..natoms-1
            index = {j: i for i, j in enumerate(atoms.ids)}

            # Convert to indices
            self.neighbors = {}
            for i, js in tmp_neighbors.items():
                self.neighbors[index[i]] = [(index[j], bo) for j, bo in js]

            # An array to keep track of the atoms visited
            self.visited = [False] * atoms.n_atoms

            # And generate the SMARTS string
            self.smarts = self._recurse(atom=0)

    def _recurse(self, atom):
        """Helper method for recursion to generate SMARTS.

        Parameters
        ----------
        atom : int
            Index of the atom to work with.
        """
        symbol = self.symbols[atom]
        self.visited[atom] = True

        # How many neighbors, and how many are hydrogens
        nX = len(self.neighbors[atom])
        nH = 0
        for neighbor, bondorder in self.neighbors[atom]:
            if self.symbols[neighbor] == "H":
                nH += 1

        # Find any non-hydrogen neighbors that have not been visited,
        # and recurse.
        sub_smarts = []
        for neighbor, bondorder in self.neighbors[atom]:
            if self.symbols[neighbor] != "H" and not self.visited[neighbor]:
                if bondorder == 1:
                    bond = ""
                elif bondorder == 2:
                    bond = "="
                elif bondorder == 3:
                    bond = "#"
                elif bondorder == 5:
                    # Aromatic
                    symbol = symbol.lower()
                    bond = ":"
                else:
                    bond = "~"
                sub_smarts.append(f"{bond}{self._recurse(neighbor)}")

        smarts = f"[{symbol}H{nH}X{nX}]"
        if len(sub_smarts) > 0:
            for sub in sub_smarts[0:-1]:
                smarts += f"({sub})"
            smarts += sub_smarts[-1]

        return smarts

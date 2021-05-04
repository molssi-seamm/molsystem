# -*- coding: utf-8 -*-

"""Functions for handling SMILES"""

import logging

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise

logger = logging.getLogger(__name__)


class SMILESMixin:
    """A mixin for handling SMILES."""

    @property
    def canonical_smiles(self):
        """Return the canonical SMILES string for this object."""
        return self.to_smiles(canonical=True)

    @property
    def smarts(self):
        """Return the SMARTS string for this object."""
        return self.to_smarts()

    @property
    def smiles(self):
        """Return the SMILES string for this object."""
        return self.to_smiles()

    def to_smiles(self, canonical=False, hydrogens=False):
        """Create the SMILES string from the system.

        Parameters
        ----------
        canonical : bool = False
            Whether to create canonical SMILES
        hydrogens : bool = False
            Whether to keep H's in the SMILES string.

        Returns
        -------
        str
            The SMILES string, or (SMILES, name) if the rname is requested
        """
        logger.info("to_smiles")

        obConversion = openbabel.OBConversion()
        if canonical:
            obConversion.SetOutFormat("can")
        else:
            obConversion.SetOutFormat("smi")

        mol = self.to_OBMol()

        if hydrogens:
            obConversion.AddOption("h")
        smiles = obConversion.WriteString(mol)

        logger.info(f"smiles = '{smiles}'")

        return smiles.strip()

    def from_smiles(self, smiles, name=None):
        """Create the system from a SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string
        name : str = None
            The name of the molecule

        Returns
        -------
        None
        """

        save = self.name

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

        molfile = obConversion.WriteString(mol)
        self.from_molfile_text(molfile)
        if name is not None:
            self.name = name
        else:
            self.name = save

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

# -*- coding: utf-8 -*-

"""Interface to RDKit."""

import logging

try:
    import rdkit
except ModuleNotFoundError:
    print(
        "Please install rdkit using conda:\n"
        "     conda install -c conda-forge rdkit"
    )
    raise

logger = logging.getLogger(__name__)


class RDKitMixin:
    """A mixin for handling RDKit via its Python interface."""

    def to_RDKMol(self):
        """Return an RDKMol object for the configuration, template, or subset."""
        rdk_mol = rdkit.Chem.rdchem.RWMol()
        for atno in self.atoms.atomic_numbers:
            _ = rdk_mol.AddAtom(rdkit.Chem.rdchem.Atom(atno))

        bond_types = {1: rdkit.Chem.rdchem.BondType.SINGLE, 2: rdkit.Chem.rdchem.BondType.DOUBLE, 3: rdkit.Chem.rdchem.BondType.TRIPLE}
        index = {j: i for i, j in enumerate(self.atoms.ids)}
        for row in self.bonds.bonds():
            rdk_mol.AddBond(index[row["i"]], index[row["j"]], bond_types[row["bondorder"]])
        
        natom = len(self.atoms.atomic_numbers)
        conf = rdkit.Chem.Conformer(natom)
        for atm_idx, xyz in self.atoms.coordinates:
            conf.SetAtomPosition(atm_idx, xyz)
        
        rdk_mol.AddConformer(conf)
        rdkit.Chem.rdmolops.SanitizeMol(rdk_mol)

        return rdk_mol

    def from_OBMol(self, rdk_mol):
        """Transform an Open Babel molecule into the current object."""
        atnos = []
        Xs = []
        Ys = []
        Zs = []
        for ob_atom in rdkit.OBMolAtomIter(rdk_mol):
            atno = ob_atom.GetAtomicNum()
            atnos.append(atno)
            Xs.append(ob_atom.x())
            Ys.append(ob_atom.y())
            Zs.append(ob_atom.z())
            logger.debug(f"atom {atno} {ob_atom.x()} {ob_atom.z()} {ob_atom.z()}")

        Is = []
        Js = []
        BondOrders = []
        for ob_bond in rdkit.OBMolBondIter(rdk_mol):
            ob_i = ob_bond.GetBeginAtom()
            ob_j = ob_bond.GetEndAtom()
            i = ob_i.GetIdx()
            j = ob_j.GetIdx()
            bondorder = ob_bond.GetBondOrder()
            Is.append(i)
            Js.append(j)
            BondOrders.append(bondorder)
            logger.debug(f"bond {i} - {j} {bondorder}")

        self.clear()
        ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)
        i = [ids[x - 1] for x in Is]
        j = [ids[x - 1] for x in Js]
        self.bonds.append(i=i, j=j, bondorder=BondOrders)

        return self

    def find_substructures(self, template):
        """Find the substructures matching the template.

        Parameters
        ----------
        template : str, _Configuration, _Template, or _Subset
            The template, which may be a SMARTS string, or a molecular object.

        Returns
        -------
        [[int]]
            Lists of atom ids for matches.
        """
        if isinstance(template, str):
            smarts = template
        else:
            smarts = self.smiles

        rdk_mol = self.to_OBMol()

        pattern = rdkit.OBSmartsPattern()
        pattern.Init(smarts)
        pattern.Match(rdk_mol)
        maplist = pattern.GetUMapList()

        return [x for x in maplist]

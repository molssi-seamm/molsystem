# -*- coding: utf-8 -*-

"""Interface to RDKit."""

import logging

try:
    import rdkit
except ModuleNotFoundError:
    print(
        "Please install rdkit using conda:\n" "     conda install -c conda-forge rdkit"
    )
    raise
import rdkit.Chem

logger = logging.getLogger(__name__)


class RDKitMixin:
    """A mixin for handling RDKit via its Python interface."""

    def to_RDKMol(self):
        """Return an RDKMol object for the configuration, template, or subset."""
        rdk_mol = rdkit.Chem.rdchem.RWMol()
        for atno in self.atoms.atomic_numbers:
            rdk_mol.AddAtom(rdkit.Chem.rdchem.Atom(atno))

        bond_types = {
            1: rdkit.Chem.rdchem.BondType.SINGLE,
            2: rdkit.Chem.rdchem.BondType.DOUBLE,
            3: rdkit.Chem.rdchem.BondType.TRIPLE,
        }
        index = {j: i for i, j in enumerate(self.atoms.ids, start=1)}
        for row in self.bonds.bonds():
            rdk_mol.AddBond(
                index[row["i"]], index[row["j"]], bond_types[row["bondorder"]]
            )

        natom = len(self.atoms.atomic_numbers)
        conf = rdkit.Chem.Conformer(natom)
        for atm_idx, xyz in enumerate(self.atoms.coordinates, start=1):
            conf.SetAtomPosition(atm_idx, xyz)

        rdk_mol.AddConformer(conf)
        rdkit.Chem.rdmolops.SanitizeMol(rdk_mol)

        return rdk_mol

    def from_RDKMol(self, rdk_mol):
        """Transform an RDKit molecule into the current object."""

        atnos = []
        for rdk_atom in rdk_mol.GetAtoms():
            atnos.append(rdk_atom.GetAtomicNum())
            logger.debug(f"atom {atnos}")

        # TODO: Generalize to handling multiple conformers
        # in a rdk_mol object, if necessary
        Xs = []
        Ys = []
        Zs = []
        for rdk_conf in rdk_mol.GetConformers():
            for atom_coords in rdk_conf.GetPositions():
                Xs.append(atom_coords[0])
                Ys.append(atom_coords[1])
                Zs.append(atom_coords[2])
                logger.debug(f"{atom_coords[0]} {atom_coords[1]} {atom_coords[2]}")

        Is = []
        Js = []
        BondOrders = []
        bond_types = {
            rdkit.Chem.rdchem.BondType.SINGLE: 1,
            rdkit.Chem.rdchem.BondType.DOUBLE: 2,
            rdkit.Chem.rdchem.BondType.TRIPLE: 3,
        }
        for rdk_bond in rdk_mol.GetBonds():
            i = rdk_bond.GetBeginAtom().GetIdx()
            j = rdk_bond.GetEndAtom().GetIdx()
            bondorder = bond_types[rdk_bond.GetBondType()]
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


#    def find_substructures(self, template):
#        """Find the substructures matching the template.
#
#        Parameters
#        ----------
#        template : str, _Configuration, _Template, or _Subset
#            The template, which may be a SMARTS string, or a molecular object.
#
#        Returns
#        -------
#        [[int]]
#            Lists of atom ids for matches.
#        """
#        if isinstance(template, str):
#            smarts = template
#        else:
#            smarts = self.smiles
#
#        rdk_mol = self.to_OBMol()
#
#        pattern = rdkit.OBSmartsPattern()
#        pattern.Init(smarts)
#        pattern.Match(rdk_mol)
#        maplist = pattern.GetUMapList()
#
#        return [x for x in maplist]

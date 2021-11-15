# -*- coding: utf-8 -*-

"""Interface to RDKit."""

import logging

try:
    from rdkit import Chem
except ModuleNotFoundError:
    print(
        "Please install rdkit using conda:\n" "     conda install -c conda-forge rdkit"
    )
    raise

logger = logging.getLogger(__name__)


class RDKitMixin:
    """A mixin for handling RDKit via its Python interface."""

    def to_RDKMol(self):
        """Return an RDKMol object for the configuration, template, or subset."""
        index = {}
        indices = []
        rdk_mol = Chem.RWMol()
        for atno, _id in zip(self.atoms.atomic_numbers, self.atoms.ids):
            idx = rdk_mol.AddAtom(Chem.Atom(atno))
            index[_id] = idx
            indices.append(idx)

        bond_types = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
        }
        for row in self.bonds.bonds():
            rdk_mol.AddBond(
                index[row["i"]], index[row["j"]], bond_types[row["bondorder"]]
            )

        natom = self.atoms.n_atoms
        conf = Chem.Conformer(natom)
        for idx, xyz in zip(indices, self.atoms.coordinates):
            conf.SetAtomPosition(idx, xyz)

        rdk_mol.AddConformer(conf)
        Chem.rdmolops.SanitizeMol(rdk_mol)

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
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
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

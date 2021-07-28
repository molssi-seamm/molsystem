# -*- coding: utf-8 -*-

"""Interface to openbabel."""

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


class OpenBabelMixin:
    """A mixin for handling OpenBabel via its Python interface."""

    def to_OBMol(self):
        """Return an OBMol object for the configuration, template, or subset."""
        ob_mol = openbabel.OBMol()
        for atno, xyz in zip(self.atoms.atomic_numbers, self.atoms.coordinates):
            ob_atom = ob_mol.NewAtom()
            ob_atom.SetAtomicNum(atno)
            ob_atom.SetVector(*xyz)

        # 1-based indices in ob.
        index = {j: i for i, j in enumerate(self.atoms.ids, start=1)}
        for row in self.bonds.bonds():
            ob_mol.AddBond(index[row["i"]], index[row["j"]], row["bondorder"])

        return ob_mol

    def from_OBMol(self, ob_mol):
        """Transform an Open Babel molecule into the current object."""
        atnos = []
        Xs = []
        Ys = []
        Zs = []
        for ob_atom in openbabel.OBMolAtomIter(ob_mol):
            atno = ob_atom.GetAtomicNum()
            atnos.append(atno)
            Xs.append(ob_atom.x())
            Ys.append(ob_atom.y())
            Zs.append(ob_atom.z())
            logger.debug(f"atom {atno} {ob_atom.x()} {ob_atom.z()} {ob_atom.z()}")

        Is = []
        Js = []
        BondOrders = []
        for ob_bond in openbabel.OBMolBondIter(ob_mol):
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

        ob_mol = self.to_OBMol()

        pattern = openbabel.OBSmartsPattern()
        pattern.Init(smarts)
        pattern.Match(ob_mol)
        maplist = pattern.GetUMapList()

        return [x for x in maplist]

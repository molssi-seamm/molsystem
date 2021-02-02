# -*- coding: utf-8 -*-

"""Interface to openbabel."""

import logging

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        'Please install openbabel using conda:\n'
        '     conda install -c conda-forge openbabel'
    )
    raise

logger = logging.getLogger(__name__)


class OpenBabelMixin:
    """A mixin for handling OpenBabel via its Python interface."""

    def to_OBMol(self):
        """Return an OBMol object for the configuration, template, or subset.
        """
        ob_mol = openbabel.OBMol()
        for atno in self.atoms.atomic_numbers:
            ob_atom = ob_mol.NewAtom()
            ob_atom.SetAtomicNum(atno)

        # 1-based indices in ob.
        index = {j: i for i, j in enumerate(self.atoms.ids, start=1)}
        for row in self.bonds.bonds():
            ob_mol.AddBond(index[row['i']], index[row['j']], row['bondorder'])

        return ob_mol

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

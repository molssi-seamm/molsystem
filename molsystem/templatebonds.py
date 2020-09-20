# -*- coding: utf-8 -*-

import logging
from typing import Any, Dict, TypeVar

from molsystem.table import _Table as Table
"""A dictionary-like object for holding bonds

Based on tables in an SQLite database.
"""

System_tp = TypeVar("System_tp", "System", "Bonds", None)

logger = logging.getLogger(__name__)


class _Templatebonds(Table):
    """The Bonds class holds arrays of attributes describing bonds

    Three attributes are required: the atoms 'i' and 'j' of the bond, and the
    bond order, which defaults to a single bond.

    Since there are two possible representations for a bond, i-j and j-i,
    the bonds are stored with i < j in order to make the indexing unique.

    In order to handle changes in bonding (and numbers of atoms) the system
    class has several tables covering the bonds. See the main documentation for
    a more detailed description. The key tables are:

    template -- a simple list of templates
    subset -- an instantiation of a template, connected with one or more
              configurations of the system.
    subset_atom -- connects the subset to the atoms, and optionally connects
                   the atom to a template atom.
    templateatom -- holds the atoms in a template, if present.
    templatebond -- holds the bonds between the template atoms.

    Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Bonds can be added ('append') or removed ('delete').
    """

    def __init__(self, system: System_tp, table: str = 'templatebond') -> None:

        super().__init__(system, table)

        self._templates = system['template']
        self._templateatoms = system['templateatom']

    def append(self, template: int = None, **kwargs: Dict[str, Any]) -> None:
        """Append one or more bonds

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used. It is an error if there is
        not a field corrresponding to a key.

        Parameters
        ----------
        template : int = None
            The template to add the bonds to.
        kwargs : Dict[str, Any]
            The attributes and their values for the bonds to append.

        Returns
        -------
        None
        """
        # Check keys and lengths of added bonds
        if 'i' not in kwargs or 'j' not in kwargs:
            raise KeyError("The atoms i & j are required!")

        if template is None:
            template = self._templates.current_template

        # How many new rows there are
        n_rows, lengths = self._get_n_rows(**kwargs)

        i = kwargs.pop('i')
        j = kwargs.pop('j')

        # Need lists for the zip below to work.
        if lengths['i'] == 0:
            i = [i] * n_rows
        if lengths['i'] == 1 and n_rows > 1:
            i = [i[0]] * n_rows
        if lengths['j'] == 0:
            j = [j] * n_rows
        if lengths['j'] == 1 and n_rows > 1:
            j = [j[0]] * n_rows

        # The list of atom ids in this template, so that we can check the atoms
        ids = self._templateatoms.atom_ids(template=template)

        i2 = []
        j2 = []
        for i_, j_ in zip(i, j):
            if not isinstance(i_, int) or not isinstance(j_, int):
                raise TypeError(
                    f"'i={i_}' and 'j={j_}', the atom indices, must be "
                    "integers"
                )
            if i_ not in ids:
                raise ValueError(f'Atom i ({i_}) is not in the template.')
            if j_ not in ids:
                raise ValueError(f'Atom j ({j_}) is not in the template.')
            # Ensure that i < j
            if i_ < j_:
                i2.append(i_)
                j2.append(j_)
            else:
                i2.append(j_)
                j2.append(i_)

        super().append(i=i2, j=j2, **kwargs)

    def bonds(self, template: int = None) -> int:
        """The an iterator over the bonds in the template

        Parameters
        ----------
        template : int = None
            The template to use, by default the current template.

        Returns
        -------
        n_bonds : int
            The number of bonds in the template
        """
        if template is None:
            template = self._templates.current_template

        return self.db.execute(
            f'SELECT * FROM {self._table}'
            '  WHERE i IN (SELECT id FROM templateatom WHERE template = ?)',
            (template,)
        )

    def n_bonds(self, template: int = None) -> int:
        """The number of bonds in the template

        Parameters
        ----------
        template : int = None
            The template to use, by default the current template.

        Returns
        -------
        n_bonds : int
            The number of bonds in the template
        """
        if template is None:
            template = self._templates.current_template

        self.cursor.execute(
            f'SELECT COUNT(*) FROM {self._table}'
            '  WHERE i IN (SELECT id FROM templateatom WHERE template = ?)'
            '    AND j IN (SELECT id FROM templateatom WHERE template = ?)',
            (template, template)
        )
        return self.cursor.fetchone()[0]

# -*- coding: utf-8 -*-

"""A dictionary-like object for holding templates


"""

import logging
from typing import TypeVar

from molsystem.table import _Table as Table

System_tp = TypeVar("System_tp", "System", None)
Templates_tp = TypeVar("Templates_tp", "_Templates", str, None)

logger = logging.getLogger(__name__)


class _Template(Table):
    """The Templates class works with SQLite tables describing templates

    See the main documentation of SEAMM for a detailed description of the
    database scheme underlying the system and hence these templates. The
    following tables handle templates:

    template -- A simple table of all of the templates, giving a name and type.

    templateatom -- An optional list of template atoms associated with a
        template.

    templatecoordinates -- An optional set of coordinates for the template
        atoms. These are in Cartesian coordinates since the templates are
        molecular, not periodic, in nature.

    templatebond -- Bonds linking template atoms i & j if present.

    Templates can be added ('append') or removed ('delete').
    """

    def __init__(self, system: System_tp, tablename: str = 'template') -> None:
        super().__init__(system, tablename)

        self._current_template = 1

    @property
    def current_template(self):
        """The template that is the default for the moment."""
        return self._current_template

    @current_template.setter
    def current_template(self, value):
        self.cursor.execute(
            f'SELECT COUNT(*) FROM {self.table} WHERE id = ?', (value,)
        )
        if self.cursor.fetchone()[0] == 0:
            raise KeyError(f"Template '{value}' does not exist.")
        self._current_template = value

    def append_to_system(
        self,
        template,
        n_copies=1,
        coordinates=None,
        configuration=None,
        create_subsets=True
    ):
        """ Append one or more copies of a template to the system.

        Parameters
        ----------

        """
        pass

    def create(self, name, type_='general', atnos=None, bonds=None):
        """Create a new template.

        Parameters
        ----------
        name : str
            The name of the template.
        type_ : str = 'general'
            The type of template, e.g. 'all', 'molecule', 'residue'
        atnos : [int] = None
            The atomic numbers of the template atoms (optional)
        bond : [(int, int, int)]
            The bonds as (i, j, order) (optional)

        Returns
        -------
        int
            The template id.
        """
        if self.exists(name, type_=type_):
            raise KeyError(f"The template '{name}' of type '{type_}' exists.")

        tid = self.append(name=name, type=type_)[0]

        if atnos is not None:
            tatom_ids = self.system.templateatoms.append(
                atno=atnos, template=tid
            )

            if bonds is not None:
                iatoms = []
                jatoms = []
                orders = []
                for i, j, order in bonds:
                    iatoms.append(tatom_ids[i])
                    jatoms.append(tatom_ids[j])
                    orders.append(order)
                self.system.templatebonds.append(
                    template=tid, i=iatoms, j=jatoms, bondorder=orders
                )
        return tid

    def exists(self, name, type_='general'):
        """Return if the template exists given the name and type.

        Parameters
        ----------
        name : str
            The name of the template.
        type_ : str = 'general'
            The type of template, e.g. 'all', 'molecule', 'residue'

        Returns
        -------
        bool
            True if it exists; False otherwise.
        """
        self.cursor.execute(
            f'SELECT COUNT(*) FROM {self.table} WHERE "name" = ? '
            'AND "type" = ?', (name, type_)
        )
        row = self.cursor.fetchone()
        return row[0] == 1

    def find(self, name, type_='general', create=False):
        """Find a single template given the name and typ.

        Parameters
        ----------
        name : str
            The name of the template.
        type_ : str = 'general'
            The type of template, e.g. 'all', 'molecule', 'residue'
        create : bool = False
            Create the template if it does not exist.

        Returns
        -------
        int
            The id of the template.

        Raises
        ------
        KeyError
           If the template does not exist and 'create' is not requested.
        """
        self.cursor.execute(
            f'SELECT id FROM {self.table} WHERE "name" = ? AND "type" = ?',
            (name, type_)
        )
        row = self.cursor.fetchone()
        if row is None:
            if create:
                return self.create(name, type_=type_)
            else:
                raise KeyError(
                    f"There is no template '{name}' of type '{type_}'."
                )
        return row[0]

    def set_current_template(self, name, type_='general'):
        """Set the current template given the name and optionally type.

        Parameters
        ----------
        name : str
            The name of the template. The name/type pair must be unique.
        type_ : str = 'general'
            The type of template.

        Returns
        -------
        None
        """
        id = self.find(name, type_=type_)
        self._current_template = id

    def templates(self, name=None, type_=None):
        """Return an itereator over the given templates.

        Parameters
        ----------
        name : str
            The name of the template.
        type_ : str = 'general'
            The type of template, e.g. 'all', 'molecule', 'residue'

        Returns
        -------
        sqlite.cursor
            The iterator over rows.
        """
        if name is None:
            if type_ is None:
                return self.db.execute(f'SELECT * FROM {self.table}')
            else:
                return self.db.execute(
                    f'SELECT * FROM {self.table} WHERE "type" = ?', (type_,)
                )
        elif type_ is None:
            return self.db.execute(
                f'SELECT * FROM {self.table} WHERE "name" = ?', (name,)
            )
        else:
            return self.db.execute(
                f'SELECT * FROM {self.table} WHERE "name" = ? AND "type" = ?',
                (name, type_)
            )

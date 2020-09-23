# -*- coding: utf-8 -*-

"""A dictionary-like object for holding atoms for the templates
"""

from itertools import zip_longest
import logging
from typing import Any, Dict, TypeVar

from molsystem.atoms import _Atoms as Atoms
from molsystem.column import _Column as Column
from molsystem.table import _Table as Table

System_tp = TypeVar("System_tp", "System", None)
Atoms_tp = TypeVar("Atoms_tp", "_Atoms", str, None)

logger = logging.getLogger(__name__)


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Templateatoms(Atoms):
    """The Atoms class holds arrays of attributes describing atoms

    This is a bit complicated due to the separation of the actual atoms and the
    coordinates, which depend on the configuration. Also, the list of atoms can
    itself be time-dependent, and is controlled by the subset 'all'.

    Atoms can be added ('append') or removed ('delete').
    """

    def __init__(
        self,
        system: System_tp,
        atom_tablename='templateatom',
        coordinates_tablename='templatecoordinates',
    ) -> None:

        self._system = system
        self._atom_tablename = atom_tablename
        self._coordinates_tablename = coordinates_tablename

        self._atom_table = Table(system, self._atom_tablename)
        self._coordinates_table = Table(system, self._coordinates_tablename)
        self._templates = system['template']

    def __getitem__(self, key) -> Any:
        """Allow [] to access the data!"""
        if key in self._atom_table.attributes:
            sql = f'WHERE template = {self.current_template}'
            return Column(self._atom_table, key, where=sql)
        elif key in self._coordinates_table.attributes:
            where = (
                "WHERE templateatom in ("
                f"     SELECT id FROM {self._atom_tablename}"
                f"      WHERE template = {self.current_template}"
                ")"
            )
            return Column(self._coordinates_table, key, where=where)
        else:
            raise KeyError(f"'{key}' not in template atoms")

    @property
    def current_template(self):
        """The current template in use."""
        return self._templates.current_template

    @current_template.setter
    def current_template(self, value):
        self._templates.current_template = value

    @property
    def n_atoms(self) -> int:
        """The number of atoms *in the current* template."""
        self.cursor.execute(
            f"SELECT COUNT(*) FROM {self._atom_tablename} WHERE template = ?",
            (self.current_template,)
        )
        return self.cursor.fetchone()[0]

    @property
    def attributes(self) -> Dict[str, Any]:
        """The definitions of the attributes.
        Combine the attributes of the atom and coordinates tables to
        make it look like a single larger table.
        """

        result = self._atom_table.attributes

        for key, value in self._coordinates_table.attributes.items():
            if key != 'templateatom':  # ignore foreign key linking tables
                result[key] = value

        return result

    @property
    def coordinate_system(self):
        """The type of coordinates: 'fractional' or 'Cartesian'"""
        return 'Cartesian'

    @coordinate_system.setter
    def coordinate_system(self, value):
        raise RuntimeError('Templates can only use Cartesian coordinates.')

    def append(self, n=None, **kwargs: Dict[str, Any]) -> None:
        """Append one or more atoms

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Need to handle the elements specially. Can give atomic numbers,
        # or symbols. By construction the references to elements are identical
        # to their atomic numbers.

        if 'symbol' in kwargs:
            symbols = kwargs.pop('symbol')
            kwargs['atno'] = self.to_atnos(symbols)

        # How many new rows there are
        n_rows, lengths = self._get_n_rows(**kwargs)

        if n is not None:
            if n_rows != 1 and n_rows != n:
                raise RuntimeError(
                    f"Requested number of template atoms ({n}) is not "
                    f"compatible with the length of the data ({n_rows})."
                )
            n_rows = n

        # Fill in the atom table
        data = {}
        for column in self._atom_table.attributes:
            if column != 'id' and column in kwargs:
                data[column] = kwargs.pop(column)
        if 'template' not in data:
            data['template'] = self.current_template

        ids = self._atom_table.append(n=n_rows, **data)

        # Now append to the coordinates table, but only if needed.
        data = {}
        for column in self._coordinates_table.attributes:
            if column != 'templateatom' and column in kwargs:
                data[column] = kwargs.pop(column)

        if len(data) > 0:
            data['templateatom'] = ids
            self._coordinates_table.append(n=n_rows, **data)

        return ids

    def atomic_numbers(self, template: int = None) -> [int]:
        """The atomic numbers of the atoms in a template.

        Parameters
        ----------
        template : int = None
            Which template, defaulting to the current template.

        Returns
        -------
        ids : [int]
            The ids of the atoms in the template.
        """
        if template is None:
            template = self.current_template
        return [
            x[0] for x in self.db.execute(
                f"SELECT atno FROM {self._atom_tablename} WHERE template = ?",
                (template,)
            )
        ]

    def atom_ids(self, template: int = None) -> [int]:
        """The ids of the atoms in a template.

        Parameters
        ----------
        template : int = None
            Which template, defaulting to the current template.

        Returns
        -------
        ids : [int]
            The ids of the atoms in the template.
        """
        if template is None:
            template = self.current_template
        return [
            x[0] for x in self.db.execute(
                f"SELECT id FROM {self._atom_tablename} WHERE template = ?",
                (template,)
            )
        ]

    def atoms(self, *args, template=None):
        """Get an iterator over atoms in the template.

        Parameters
        ----------
        args : str, int or float
            SQL restrictions for a WHERE statement, each argument being
            one word, e.g. "atno" "=" 5
        template : int = None
            Which template, defaulting to the current template.

        Returns
        -------
        SQLite3.Cursor :
            The cursor containing the result.
        """
        atom_tbl = self._atom_tablename
        coord_tbl = self._coordinates_tablename
        atom_columns = [*self._atom_table.attributes]
        coord_columns = [*self._coordinates_table.attributes]
        coord_columns.remove('templateatom')

        columns = [f'{atom_tbl}.{x}' for x in atom_columns]
        columns += [f'{coord_tbl}.{x}' for x in coord_columns]
        column_defs = ', '.join(columns)

        if template is None:
            template = self.current_template

        sql = (
            f'SELECT {column_defs} FROM {atom_tbl}, {coord_tbl}'
            f' WHERE {coord_tbl}.templateatom = {atom_tbl}.id'
            f'   AND {atom_tbl}.template = ?'
        )
        if len(args) == 0:
            return self.db.execute(sql, (template,))

        parameters = [template]
        for col, op, value in grouped(args, 3):
            if op == '==':
                op = '='
            sql += f' AND "{col}" {op} ?'
            parameters.append(value)

        return self.db.execute(sql, parameters)

    def remove(self, *args, template=None):
        """Remove atoms in the template.

        Parameters
        ----------
        args : str, int or float
            SQL restrictions for a WHERE statement, each argument being
            one word, e.g. "atno" "=" 5
        template : int = None
            Which template, defaulting to the current template.

        Returns
        -------
        None
        """
        if template is None:
            template = self.current_template

        sql = f'DELETE FROM {self.table} WHERE template = ?'

        parameters = [template]
        for col, op, value in grouped(args, 3):
            if op == '==':
                op = '='
            sql += f' AND "{col}" {op} ?'
            parameters.append(value)

        self.db.execute(sql, parameters)
        self.db.commit()

# -*- coding: utf-8 -*-

"""A dictionary-like object for holding atoms

For efficiency the atom-data is stored as arrays -- numpy if possible, python
lists otherwise, along with metadata to define the attributes ("columns").

In some ways this a bit like pandas; however, we also need more control than a
simple pandas dataframe provides....
"""

import collections.abc
from itertools import zip_longest
import logging
from typing import Any, Dict, TypeVar

import numpy
import pandas

from molsystem.column import _Column as Column
from molsystem.table import _Table as Table

System_tp = TypeVar("System_tp", "System", None)
Atoms_tp = TypeVar("Atoms_tp", "_Atoms", str, None)

logger = logging.getLogger(__name__)


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Atoms(collections.abc.MutableMapping):
    """The Atoms class holds arrays of attributes describing atoms.

    In order to handle changes in bonding (and numbers of atoms) the
    system class has several tables covering the atoms and bonds. See
    the main documentation for a more detailed description. The key
    tables are:

    configuration -- a list of configurations in e.g. a trajectory
    configuration_subset -- connnects configurations with subset

    template -- a simple list of templates
    templateatom -- holds the atoms in a template, if present.
    templatebond -- holds the bonds between the template atoms, if present

    subset -- an instantiation of a template, connected with one or more
              configurations of the system.
    subset_atom -- connects the subset to the atoms, and optionally connects
                   the atom to a template atom.

    atom -- the nonvarying part of the description of the atoms
    coordinates -- the varying part of the description of atoms

    The description of the atoms is split into two parts: the identity
    of the atom in the 'atom' table; and the coordinates and other
    varying attributes in the table 'coordinates' as indicated above.

    This class handles the situation where the number and type of atoms
    changes from configuration to configuration, as it would for example
    in a grand canonical simulation or simulating surface deposition.

    The configuration if either given explicitly or the "current configuration"
    is used. The current configuration is set and stored in the system object.

    Given this configuration, this class provides the illusion that
    there is a single set of atoms in a single table containing both the
    identity of the atoms and their coordinates.

    Underneath, this is handled by a subset "all", which is a special
    subset that contains all of the atoms in the system. This subset is
    defined by one or more templates "all" which are instantiated in the
    subset table and connected to configurations via the
    configuration_subset table.

    If the number of atoms in the system is fixed and there are no bonds
    or they don't change (the case in classical MD), there will be one
    template "all" instantiated as one subset. All configurations will
    point to this one subset. If there are no bonds in the system, then
    this is all that is needed. However if there are bonds, the atoms
    involved in bonds need to be in the templateatom table and the bonds
    entered in the templatebond table.

    If the bonds change over time, then for each bond configuration
    there needs to be a distinct "all" template, with the appropriate
    templateatoms and templatebonds. This template is instantiated as a
    new subset and connected to one or more configurations.

    If the number of atoms changes over time, each new set of atoms also
    requires a new template and subset connected to the atoms in that
    configuration. A subset may persist over a number of configurations,
    but each time the atoms change a new subset is required.

    If both the atoms and the bonds change, then a corresponding subset
    is needed each time either changes.
    """

    def __init__(
        self,
        system: System_tp,
        atom_tablename='atom',
        coordinates_tablename='coordinates',
    ) -> None:

        self._system = system
        self._atom_tablename = atom_tablename
        self._coordinates_tablename = coordinates_tablename

        self._atom_table = Table(system, self._atom_tablename)
        self._coordinates_table = Table(system, self._coordinates_tablename)

    def __enter__(self) -> Any:
        self.system.__enter__()
        return self

    def __exit__(self, etype, value, traceback) -> None:
        self.system.__exit__(etype, value, traceback)

    def __getitem__(self, key) -> Any:
        """Allow [] to access the data!"""
        return self.get_column(key)

    def __setitem__(self, key, value) -> None:
        """Allow x[key] access to the data"""
        column = self[key]
        column[0:] = value

    def __delitem__(self, key) -> None:
        """Allow deletion of keys"""
        if key in self._atom_table.attributes:
            del self._atom_table[key]
        elif key in self._coordinates_table.attributes:
            del self._coordinates_table[key]

    def __iter__(self) -> iter:
        """Allow iteration over the object"""
        return iter([*self.attributes.keys()])

    def __len__(self) -> int:
        """The len() command"""
        return len(self.attributes)

    def __repr__(self) -> str:
        """The string representation of this object"""
        df = self.to_dataframe()
        return repr(df)

    def __str__(self) -> str:
        """The pretty string representation of this object"""
        df = self.to_dataframe()
        return str(df)

    def __contains__(self, item) -> bool:
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.

        tmp_item = item.strip('"')
        return tmp_item in self.attributes

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        if self._atom_table != other._atom_table:
            return False
        if self._coordinates_table != other._coordinates_table:
            return False
        return True

    @property
    def system(self):
        """The system that we belong to."""
        return self._system

    @property
    def db(self):
        """The database connection."""
        return self.system.db

    @property
    def cursor(self):
        """The database connection."""
        return self.system.cursor

    @property
    def table(self) -> str:
        """The name of this table."""
        return '"' + self._atom_tablename + '"'

    @property
    def current_configuration(self):
        """The configuration of the system being used"""
        return self.system.current_configuration

    @current_configuration.setter
    def current_configuration(self, value):
        self.system.current_configuration = value

    @property
    def n_rows(self) -> int:
        """The number of rows in the table."""
        self.cursor.execute(f'SELECT COUNT(*) FROM {self.table}')
        result = self.cursor.fetchone()[0]
        return result

    @property
    def attributes(self) -> Dict[str, Any]:
        """The definitions of the attributes.
        Combine the attributes of the atom and coordinates tables to
        make it look like a single larger table.
        """

        result = self._atom_table.attributes

        for key, value in self._coordinates_table.attributes.items():
            if key != 'atom':  # atom key links the tables together, so ignore
                result[key] = value

        return result

    @property
    def coordinate_system(self):
        """The type of coordinates: 'fractional' or 'Cartesian'"""
        return self._system.coordinate_system

    @coordinate_system.setter
    def coordinate_system(self, value):
        self._system.coordinate_system = value

    @property
    def version(self):
        return self.system.version

    def add_attribute(
        self,
        name: str,
        coltype: str = 'float',
        default: Any = None,
        notnull: bool = False,
        index: bool = False,
        pk: bool = False,
        references: str = None,
        values: Any = None,
        configuration_dependent: bool = False
    ) -> None:
        """Adds a new attribute.

        If the default value is None, you must always provide values wherever
        needed, for example when adding a row.

        Parameters
        ----------
            name : str
                the name of the attribute.
            coltype : str = 'float'
                the type of the attribute (column). Must be one of 'int',
                'float', 'str' or 'byte'. Defaults to 'float'.
            default : Any
                the default value for the attribute if no value is given.
            notnull : bool = False
                whether the value must be non-null
            index :  bool = False
                whether to create an index on the column
            pk : bool
                whether the column is the primay key
            references : str = None
                If the column reference another table, i.e. is a FK
            values : Any
                either a single value or a list of values length 'nrows' of
                values to fill the column.
            configuration_dependent : bool = False
                whether the attribute belongs with the coordinates (True)
                or atoms (False)

        Returns
        -------
            None
        """

        if configuration_dependent:
            self._coordinates_table.add_attribute(
                name, coltype, default, notnull, index, pk, references, values
            )
        else:
            self._atom_table.add_attribute(
                name, coltype, default, notnull, index, pk, references, values
            )

    def append(self, configuration=None, **kwargs: Dict[str, Any]) -> None:
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

        # Fill in the atom table
        data = {}
        for column in self._atom_table.attributes:
            if column != 'id' and column in kwargs:
                data[column] = kwargs.pop(column)

        if len(data) == 0:
            data['atno'] = [None] * n_rows

        ids = self._atom_table.append(n=n_rows, **data)

        # Now append to the coordinates table
        if configuration is None:
            configuration = self.current_configuration
        data = {'configuration': configuration, 'atom': ids}
        for column in self._coordinates_table.attributes:
            if column != 'id' and column in kwargs:
                data[column] = kwargs.pop(column)

        self._coordinates_table.append(n=n_rows, **data)

        # And to the subset 'all'
        subset_atom = self.system['subset_atom']
        subset_atom.append(
            subset=self.system.all_subset(configuration), atom=ids
        )

        return ids

    def atoms(self, *args, configuration=None):
        """Return an iterator over the atoms."""
        atom_tbl = self._atom_tablename
        coord_tbl = self._coordinates_tablename
        atom_columns = [*self._atom_table.attributes]
        coord_columns = [*self._coordinates_table.attributes]
        coord_columns.remove('atom')

        columns = [f'{atom_tbl}.{x}' for x in atom_columns]
        columns += [f'{coord_tbl}.{x}' for x in coord_columns]
        column_defs = ', '.join(columns)

        sql = (
            f'SELECT {column_defs} FROM {atom_tbl}, {coord_tbl}'
            f' WHERE {atom_tbl}.id IN ('
            f'   SELECT atom FROM subset_atom WHERE subset = ?'
            f') AND {coord_tbl}.atom = {atom_tbl}.id'
        )
        all_subset = self.system.all_subset(configuration)
        if len(args) == 0:
            return self.db.execute(sql, (all_subset,))

        parameters = [all_subset]
        for col, op, value in grouped(args, 3):
            if op == '==':
                op = '='
            sql += f' AND "{col}" {op} ?'
            parameters.append(value)

        return self.db.execute(sql, parameters)

    def atom_ids(self, configuration=None) -> [int]:
        """The ids of the atoms the configuration."""
        return [
            x[0] for x in self.db.execute(
                "SELECT atom FROM subset_atom WHERE subset = ?",
                (self.system.all_subset(configuration),)
            )
        ]

    def atomic_numbers(self, configuration: int = None) -> [int]:
        """The atomic numbers of the atoms in the configuration."""
        return [*self['atno']]

    def coordinates(
        self,
        configuration=None,
        fractionals=True,
        in_cell=False,
        as_array=False
    ):
        """Return the coordinates optionally translated back into the principal
        unit cell.

        Parameters
        ----------
        configuration : int = None
            The configuration of interest.
        frationals : bool = True
            Return the coordinates as fractional coordinates for periodic
            systems. Non-periodic systems always use Cartesian coordinates.
        in_cell : bool, str = False
            Whether to translate the atoms into the unit cell, and if so
            whether to do so by molecule or just atoms.
        as_array : bool = False
            Whether to return the results as a numpy array or as a list of
            lists (the default).

        Returns
        -------
        abc : [N][float*3]
            The coordinates, either Cartesian or fractional
        """
        if configuration is None:
            configuration = self.current_configuration

        xyz = []
        for row in self.atoms(configuration=configuration):
            xyz.append([row['x'], row['y'], row['z']])

        periodicity = self.system.periodicity
        if periodicity == 0:
            if as_array:
                return numpy.array(xyz)
            else:
                return xyz

        cell = self.system['cell'].cell(configuration)

        if 'molecule' in in_cell:
            # Need fractionals...
            if self.system.coordinate_system == 'Cartesian':
                UVW = cell.to_fractionals(xyz, as_array=True)
            elif not isinstance(xyz, numpy.ndarray):
                UVW = numpy.array(xyz)
            else:
                UVW = xyz

            molecules = self.system.find_molecules(configuration=configuration)

            for indices in molecules:
                indices = numpy.array([i - 1 for i in indices])
                uvw_mol = numpy.take(UVW, indices, axis=0)
                center = numpy.average(uvw_mol, axis=0)
                delta = numpy.floor(center)
                uvw_mol -= delta
                numpy.put_along_axis(
                    UVW, numpy.expand_dims(indices, axis=1), uvw_mol, axis=0
                )
            if fractionals:
                if as_array:
                    return UVW
                else:
                    return UVW.tolist()
            else:
                return cell.to_cartesians(UVW, as_array=as_array)
        elif in_cell:
            # Need fractionals...
            if self.system.coordinate_system == 'Cartesian':
                UVW = cell.to_fractionals(xyz, as_array=True)
            elif not isinstance(xyz, numpy.ndarray):
                UVW = numpy.array(xyz)
            else:
                UVW = xyz
            delta = numpy.floor(UVW)
            UVW -= delta
            if fractionals:
                if as_array:
                    return UVW
                else:
                    return UVW.tolist()
            else:
                return cell.to_cartesians(UVW, as_array=as_array)
        else:
            if fractionals:
                if self.system.coordinate_system == 'Cartesian':
                    return cell.to_fractionals(xyz, as_array=as_array)
                elif as_array:
                    return numpy.array(xyz)
                else:
                    return xyz
            else:
                if self.system.coordinate_system == 'fractional':
                    return cell.to_cartesians(xyz, as_array=as_array)
                elif as_array:
                    return numpy.array(xyz)
                else:
                    return xyz

    def get_column(self, key, configuration=None) -> Any:
        """Allow [] to access the data!"""
        if configuration is None:
            configuration = self.system.current_configuration
        if key in self._atom_table.attributes:
            sql = (
                'WHERE id IN (SELECT atom FROM subset_atom '
                f'WHERE subset = {self.system.all_subset(configuration)})'
            )
            return Column(self._atom_table, key, where=sql)
        elif key in self._coordinates_table.attributes:
            where = f"WHERE configuration = {configuration}"
            return Column(self._coordinates_table, key, where=where)
        else:
            raise KeyError(f"'{key}' not in atoms")

    def to_atnos(self, symbols):
        """Convert element symbols to atomic numbers."""
        return self._system.to_atnos(symbols)

    def _get_n_rows(self, **kwargs):
        """Get the total number of rows represented in the arguments."""
        n_rows = None

        lengths = {}
        for key, value in kwargs.items():
            if key not in self:
                raise KeyError(f'"{key}" is not an attribute of the atoms.')
            length = self.length_of_values(value)
            lengths[key] = length

            if n_rows is None:
                n_rows = 1 if length == 0 else length

            if length > 1 and length != n_rows:
                if n_rows == 1:
                    n_rows = length
                else:
                    raise IndexError(
                        'key "{}" has the wrong number of values, '
                        .format(key) +
                        '{}. Should be 1 or the number of atoms ({}).'
                        .format(length, n_rows)
                    )
        return n_rows, lengths

    def length_of_values(self, values: Any) -> int:
        """Return the length of the values argument.

        Parameters
        ----------
        values : Any
            The values, which might be a string, single value, list, tuple,
            etc.

        Returns
        -------
        length : int
            The length of the values. 0 indicates a scalar
        """
        if isinstance(values, str):
            return 0
        else:
            try:
                return len(values)
            except TypeError:
                return 0

    def n_atoms(self, configuration=None) -> int:
        """The number of atoms in a configuration."""
        self.cursor.execute(
            "SELECT COUNT(*) FROM subset_atom WHERE subset = ?",
            (self.system.all_subset(configuration),)
        )
        return self.cursor.fetchone()[0]

    def symbols(self, configuration: int = None) -> [str]:
        """Convert element symbols to atomic numbers."""
        return self._system.to_symbols(self.atomic_numbers(configuration))

    def to_dataframe(self):
        """Return the contents of the table as a Pandas Dataframe."""
        data = {}
        rows = self.atoms()
        for row in rows:
            data[row[0]] = row[1:]

        columns = [x[0] for x in rows.description[1:]]
        df = pandas.DataFrame.from_dict(data, orient='index', columns=columns)

        return df


if __name__ == '__main__':  # pragma: no cover
    import numpy.random as nprand
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False) -> None:
        system = None
        atoms = _Atoms(system)
        x = nprand.uniform(low=0, high=100, size=nper)
        y = nprand.uniform(low=0, high=100, size=nper)
        z = nprand.uniform(low=0, high=100, size=nper)
        atno = numpy.array(nper * [6])
        with atoms as tmp:
            for i in range(0, nrepeat):
                tmp.append(x=x, y=y, z=z, atno=atno)
                x += 10.0
                y += 5.0
                z -= 3.0

    nrepeat = 1
    nper = 10
    t = timeit.timeit(
        "run(nper={}, nrepeat={})".format(nper, nrepeat),
        setup="from __main__ import run",
        timer=time.time,
        number=1
    )

    print("Creating {} atoms took {:.3f} s".format(nper * nrepeat, t))

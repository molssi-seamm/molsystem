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
        on_delete: str = 'cascade',
        on_update: str = 'cascade',
        values: Any = None,
        configuration_dependent: bool = False
    ) -> None:
        """Adds a new attribute.

        If the default value is None, you must always provide values wherever
        needed, for example when adding a row.

        Parameters
        ----------
        name : str
            The name of the attribute.
        coltype : str
            The type of the attribute (column). Must be one of 'int',
            'float', 'str' or 'byte'
        default : int, float, str or byte
            The default value for the attribute if no value is given.
        notnull : bool = False
            Whether the value must be non-null
        index : bool = False
            Whether to create an index on the column
        pk : bool = False
            Whether the column is the primry keys
        references : str = None
            If not null, the column is a foreign key for this table.
        on_delete : str = 'cascade'
            How to handle deletions of a foregin keys
        on_update : str = 'cascade'
            How to handle updates of a foregin key
        values : Any
            Either a single value or a list of values length 'nrows' of
            values to fill the column.
        configuration_dependent : bool = False
            Whether the attribute belongs with the coordinates (True)
            or atoms (False)

        Returns
        -------
            None
        """
        if configuration_dependent:
            self._coordinates_table.add_attribute(
                name,
                coltype=coltype,
                default=default,
                notnull=notnull,
                index=index,
                pk=pk,
                references=references,
                on_delete=on_delete,
                on_update=on_update,
                values=values
            )
        else:
            self._atom_table.add_attribute(
                name,
                coltype=coltype,
                default=default,
                notnull=notnull,
                index=index,
                pk=pk,
                references=references,
                on_delete=on_delete,
                on_update=on_update,
                values=values
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
        if isinstance(configuration, int):
            config = configuration
        else:
            config = configuration[0]

        subset_atom.append(subset=self.system.all_subset(config), atom=ids)

        return ids

    def atoms(
        self, *args, subset=None, configuration=None, template_order=False
    ):
        """Return an iterator over the atoms.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            return rows in the order of the template.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the atoms.
        """
        if subset is None:
            subset = self.system.all_subset(configuration)

        # If we are asked to use template order, see if all the atoms
        # have associated template atoms.
        if template_order:
            self.cursor.execute(
                'SELECT COUNT(*) FROM subset_atom WHERE subset = ?'
                '    AND templateatom IS NULL', (subset,)
            )
            n = self.cursor.fetchone()[0]
            if n > 0:
                raise RuntimeError(
                    'Not all of the atoms are defined in the template for '
                    f'subset {subset} - {n} are not'
                )

        atom_tbl = self._atom_tablename
        coord_tbl = self._coordinates_tablename
        atom_columns = [*self._atom_table.attributes]
        coord_columns = [*self._coordinates_table.attributes]
        coord_columns.remove('atom')

        columns = [f'at.{x}' for x in atom_columns]
        columns += [f'co.{x}' for x in coord_columns]
        column_defs = ', '.join(columns)

        sql = (
            f'SELECT {column_defs}'
            f'  FROM {atom_tbl} as at, {coord_tbl} as co, subset_atom as sa'
            '  WHERE at.id == sa.atom AND sa.subset = ? AND co.atom = at.id'
        )

        if len(args) == 0:
            if template_order:
                sql += " ORDER BY templateatom"
            return self.db.execute(sql, (subset,))

        parameters = [subset]
        for col, op, value in grouped(args, 3):
            if op == '==':
                op = '='
            sql += f' AND "{col}" {op} ?'
            parameters.append(value)
        if template_order:
            sql += " ORDER BY templateatom"

        return self.db.execute(sql, parameters)

    def atom_ids(self, subset=None, configuration=None, template_order=False):
        """The ids of the atoms the subset or configuration.

        Parameters
        ----------
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            return rows in the order of the template.

        Returns
        -------
        [int]
            The ids of the requested atoms.
        """
        if subset is None:
            subset = self.system.all_subset(configuration)
        if template_order:
            return [
                x[0] for x in self.db.execute(
                    "SELECT atom, templateatom FROM subset_atom "
                    "WHERE subset = ? ORDER BY templateatom", (subset,)
                )
            ]
        else:
            return [
                x[0] for x in self.db.execute(
                    "SELECT atom FROM subset_atom WHERE subset = ?", (subset,)
                )
            ]

    def atomic_numbers(
        self,
        subset: int = None,
        configuration: int = None,
        template_order: bool = False
    ) -> [float]:
        """The atomic numbers of the atoms in the subset or configuration.

        Parameters
        ----------
        subset : int = None
            Get the values for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            return rows in the order of the template.

        Returns
        -------
        [int]
            The atomic numbers.
        """

        column = self.get_column(
            'atno',
            subset=subset,
            configuration=configuration,
            template_order=template_order
        )
        return [*column]

    def atomic_masses(
        self,
        subset: int = None,
        configuration: int = None,
        template_order: bool = False
    ) -> [float]:
        """The atomic masses of the atoms in the subset or configuration.

        Parameters
        ----------
        subset : int = None
            Get the values for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            return rows in the order of the template.

        Returns
        -------
        [int]
            The atomic numbers.
        """

        if 'mass' in self:
            column = self.get_column(
                'mass',
                subset=subset,
                configuration=configuration,
                template_order=template_order
            )
            return [*column]
        else:
            atnos = self.atomic_numbers(
                subset=subset,
                configuration=configuration,
                template_order=template_order
            )
            return self._system.default_masses(atnos=atnos)

    def coordinates(
        self,
        subset=None,
        configuration=None,
        template_order=False,
        fractionals=True,
        in_cell=False,
        as_array=False
    ):
        """Return the coordinates optionally translated back into the principal
        unit cell.

        Parameters
        ----------
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            return rows in the order of the template.
        fractionals : bool = True
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
        for row in self.atoms(
            subset=subset,
            configuration=configuration,
            template_order=template_order
        ):
            xyz.append([row['x'], row['y'], row['z']])

        periodicity = self.system.periodicity
        if periodicity == 0:
            if as_array:
                return numpy.array(xyz)
            else:
                return xyz

        cell = self.system['cell'].cell(configuration)

        if isinstance(in_cell, str) and 'molecule' in in_cell:
            # Need fractionals...
            if self.system.coordinate_system == 'Cartesian':
                UVW = cell.to_fractionals(xyz, as_array=True)
            elif not isinstance(xyz, numpy.ndarray):
                UVW = numpy.array(xyz)
            else:
                UVW = xyz

            molecules = self.system.find_molecules(
                configuration=configuration, as_indices=True
            )

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

    def set_coordinates(
        self,
        xyz,
        subset=None,
        configuration=None,
        template_order=False,
        fractionals=True
    ):
        """Set the coordinates to new values.

        Parameters
        ----------
        subset : int = None
            Set the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.
        template_order : bool = False
            If True, and there are template atoms associated with the atoms,
            the coordinates are in the order of the template.
        fractionals : bool = True
            The coordinates are fractional coordinates for periodic
            systems. Ignored for non-periodic systems.

        Returns
        -------
        None
        """
        if configuration is None:
            configuration = self.current_configuration

        as_array = isinstance(xyz, numpy.ndarray)

        x_column = self.get_column(
            'x',
            subset=subset,
            configuration=configuration,
            template_order=template_order
        )
        y_column = self.get_column(
            'y',
            subset=subset,
            configuration=configuration,
            template_order=template_order
        )
        z_column = self.get_column(
            'z',
            subset=subset,
            configuration=configuration,
            template_order=template_order
        )

        xs = []
        ys = []
        zs = []

        periodicity = self.system.periodicity
        coord_system = self.system.coordinate_system
        if (
            periodicity == 0 or
            (coord_system == 'Cartesian' and not fractionals) or
            (coord_system == 'fractional' and fractionals)
        ):
            if as_array:
                for x, y, z in xyz.tolist():
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            else:
                for x, y, z in xyz:
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
        else:
            cell = self.system['cell'].cell(configuration)
            if coord_system == 'fractional':
                # Convert coordinates to fractionals
                for x, y, z in cell.to_fractionals(xyz):
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            else:
                for x, y, z in cell.to_cartesians(xyz):
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
        x_column[0:] = xs
        y_column[0:] = ys
        z_column[0:] = zs

    def get_column(
        self,
        key: str,
        subset: int = None,
        configuration: int = None,
        template_order: bool = False
    ) -> Any:
        """Get a Column object with the requested data

        Parameters
        ----------
        key : str
            The attribute to get.
        subset : int = None
            Get the values for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        Column
            A Column object containing the data.
        """
        if subset is None:
            subset = self.system.all_subset(configuration)

        if key in self._atom_table.attributes:
            sql = (
                f'SELECT at.rowid, at.{key}, sa.templateatom'
                f'  FROM {self._atom_tablename} as at, subset_atom as sa'
                f' WHERE at.id = sa.atom AND sa.subset = {subset}'
            )
            if template_order:
                sql += " ORDER BY sa.templateatom"
            return Column(self._atom_table, key, sql=sql)
        elif key in self._coordinates_table.attributes:
            sql = (
                f'SELECT co.rowid, co.{key}, sa.templateatom'
                f'  FROM {self._atom_tablename} as at,'
                f'       {self._coordinates_tablename} as co,'
                '        subset_atom as sa'
                f' WHERE co.atom = at.id AND at.id = sa.atom'
                f'   AND sa.subset = {subset}'
            )
            if template_order:
                sql += " ORDER BY sa.templateatom"
            return Column(self._coordinates_table, key, sql=sql)
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

    def n_atoms(self, subset=None, configuration=None) -> int:
        """The number of atoms in a subset or configuration

        Parameters
        ----------
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        int
            Number of atoms
        """
        if subset is None:
            subset = self.system.all_subset(configuration)

        self.cursor.execute(
            "SELECT COUNT(*) FROM subset_atom WHERE subset = ?", (subset,)
        )
        return self.cursor.fetchone()[0]

    def remove(self, atoms=None, subset=None, configuration=None) -> int:
        """Delete the atoms listed, or in a subset or configuration

        Parameters
        ----------
        atoms : [int] = None
            The list of atoms to delete
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        None
        """
        if atoms is not None:
            # Delete the listed atoms and coordinates
            subset = self.system.all_subset(configuration)

            # Need to handle bonds first
            self.system.bonds.remove(atoms=atoms, subset=subset)

            parameters = [(i,) for i in atoms]
            # Coordinates
            self.db.executemany(
                "DELETE FROM coordinates WHERE atom = ?", parameters
            )

            # Atoms
            self.db.executemany("DELETE FROM atom WHERE id = ?", parameters)

            # Subset-Atoms
            self.db.executemany(
                "DELETE FROM subset_atom WHERE atom = ?", parameters
            )

            return
        if subset is None:
            subset = self.system.all_subset(configuration)
            # Bonds only if removing all atoms, i.e. subset all
            self.system.bonds.remove(subset=subset)

        # Coordinates
        self.db.execute(
            "DELETE FROM coordinates"
            " WHERE atom IN ("
            "   SELECT atom FROM subset_atom WHERE subset = ?"
            ")", (subset,)
        )

        # Atoms
        self.db.execute(
            "DELETE FROM atom"
            " WHERE id IN ("
            "   SELECT atom FROM subset_atom WHERE subset = ?"
            ")", (subset,)
        )

        # Subset-Atoms
        self.db.execute("DELETE FROM subset_atom WHERE subset = ?", (subset,))

    def symbols(self, subset=None, configuration: int = None) -> [str]:
        """The element symbols for the atoms in a subset or configuration

        Parameters
        ----------
        subset : int = None
            Get the values for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        [str]
            The element symbols
        """

        return self._system.to_symbols(
            self.atomic_numbers(subset=subset, configuration=configuration)
        )

    def to_dataframe(self):
        """Return the contents of the table as a Pandas Dataframe."""
        data = {}
        rows = self.atoms()
        for row in rows:
            data[row[0]] = row[1:]

        columns = [x[0] for x in rows.description[1:]]
        df = pandas.DataFrame.from_dict(data, orient='index', columns=columns)

        return df

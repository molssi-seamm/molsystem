# -*- coding: utf-8 -*-

"""A dictionary-like object for holding atoms.
"""

from itertools import zip_longest
import logging
from typing import Any, Dict, TypeVar

import numpy as np
import pandas

from molsystem import elements
from .column import _Column
from .table import _Table

Atoms_tp = TypeVar("Atoms_tp", "_Atoms", str, None)

logger = logging.getLogger(__name__)
# logger.setLevel("DEBUG")


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Atoms(_Table):
    """The atoms in a configuration of a system.

    Parameters
    ----------
    configuration : _Configuration
        The configuration of interest.
    """

    def __init__(self, configuration) -> None:
        self._configuration = configuration

        self._system_db = self._configuration.system_db
        self._system = None

        self._atomset = self._configuration.atomset

        self._atom_table = _Table(self.system_db, "atom")
        self._coordinates_table = _Table(self.system_db, "coordinates")
        self._velocities_table = _Table(self.system_db, "velocities")

        super().__init__(self._system_db, "atom")

    def __enter__(self) -> Any:
        """Copy the tables to a backup for a 'with' statement."""
        self.system_db["atomset_atom"].__enter__()
        self.system_db["atom"].__enter__()
        return self

    def __exit__(self, etype, value, traceback) -> None:
        """Handle returning from a 'with' statement."""
        if etype is None:
            self.configuration.version = self.configuration.version + 1
        self.system_db["atomset_atom"].__exit__(etype, value, traceback)
        return self.system_db["atom"].__exit__(etype, value, traceback)

    def __delitem__(self, key) -> None:
        """Allow deletion of keys"""
        if key in self._atom_table.attributes:
            del self._atom_table[key]
        elif key in self._coordinates_table.attributes:
            del self._coordinates_table[key]
        elif key in self._velocities_table.attributes:
            del self._velocities_table[key]

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
        return df.to_string()

    def __contains__(self, item) -> bool:
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.

        tmp_item = item.strip('"')
        return tmp_item in self.attributes

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        diffs = self.diff(other)
        return len(diffs) == 0

    @property
    def ids(self):
        """The ids of the atoms."""
        return self.get_ids()

    @property
    def asymmetric_atomic_numbers(self):
        """The atomic numbers of the asymmetric atoms.

        Returns
        -------
        [int]
            The atomic numbers.
        """
        return self.get_column_data("atno")

    @property
    def atomic_numbers(self):
        """The atomic numbers of the atoms.

        Returns
        -------
        [int]
            The atomic numbers.
        """
        if self.configuration.symmetry.n_symops == 1:
            return self.asymmetric_atomic_numbers
        else:
            result = []
            atno = self.asymmetric_atomic_numbers
            for i in self.configuration.atom_to_asymmetric_atom:
                result.append(atno[i])
            return result

    @property
    def asymmetric_atomic_masses(self):
        """The atomic masses of the atoms.

        Returns
        -------
        [int]
            The atomic numbers.
        """

        if "mass" in self:
            result = self.get_column_data("mass")
        else:
            atnos = self.atomic_numbers
            result = elements.masses(atnos)
        return result

    @property
    def atomic_masses(self):
        """The atomic masses of the atoms.

        Returns
        -------
        [int]
            The atomic numbers.
        """
        if self.configuration.symmetry.n_symops == 1:
            return self.asymmetric_atomic_masses
        else:
            result = []
            mass = self.asymmetric_atomic_masses
            for i in self.configuration.atom_to_asymmetric_atom:
                result.append(mass[i])
            return result

    @property
    def atomset(self):
        """The atomset for these atoms."""
        return self._atomset

    @property
    def atom_generators(self):
        """The symmetry operations that create the symmetric atoms."""
        return self.configuration.symmetry.atom_generators

    @property
    def attributes(self) -> Dict[str, Any]:
        """The definitions of the attributes.
        Combine the attributes of the atom, coordinates, and velocities tables to
        make it look like a single larger table.
        """

        result = self._atom_table.attributes

        for key, value in self._coordinates_table.attributes.items():
            if key != "atom":  # atom key links the tables together, so ignore
                result[key] = value

        for key, value in self._velocities_table.attributes.items():
            if key != "atom":  # atom key links the tables together, so ignore
                result[key] = value

        return result

    @property
    def configuration(self):
        """Return the configuration."""
        return self._configuration

    @property
    def coordinates(self):
        """The coordinates as list of lists."""
        return self.get_coordinates()

    @coordinates.setter
    def coordinates(self, xyz):
        """The coordinates as list of lists."""
        return self.set_coordinates(xyz)

    @property
    def cursor(self):
        """The database connection."""
        return self.system_db.cursor

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def group(self):
        """The space or point group of the configuration."""
        return self.configuration.symmetry.group

    @group.setter
    def group(self, value):
        self.configuration.symmetry.group = value

    @property
    def have_velocities(self):
        """Whether there are any velocities for this configuration."""
        sql = (
            "SELECT COUNT(*)"
            " FROM atomset_atom AS aa, velocities AS ve "
            "WHERE aa.atomset = ?"
            "  AND ve.atom = aa.atom AND ve.configuration = ?"
        )
        parameters = [self.atomset, self.configuration.id]
        self.cursor.execute(sql, parameters)
        return self.cursor.fetchone()[0] > 0

    @property
    def loglevel(self):
        """The logging level for this module."""
        result = logger.getEffectiveLevel()
        tmp = logging.getLevelName(result)
        if "Level" not in tmp:
            result = tmp
        return result

    @loglevel.setter
    def loglevel(self, value):
        logger.setLevel(value)

    @property
    def names(self):
        """The names of the atoms."""
        return self.get_names()

    @property
    def n_asymmetric_atoms(self) -> int:
        """The number of symmetry-unique atoms in this configuration."""
        self.cursor.execute(
            "SELECT COUNT(*) FROM atomset_atom WHERE atomset = ?", (self.atomset,)
        )
        result = self.cursor.fetchone()[0]
        return result

    @property
    def n_atoms(self) -> int:
        """The number of atoms in this configuration."""
        if self.configuration.symmetry.n_symops == 1:
            return self.n_asymmetric_atoms
        else:
            n_atoms = 0
            for symops in self.atom_generators:
                n_atoms += len(symops)
            return n_atoms

    @property
    def n_symops(self):
        """The number of symmetry operations in the group."""
        return self.configuration.symmetry.n_symops

    @property
    def asymmetric_symbols(self):
        """The element symbols for the atoms in this configuration.

        Returns
        -------
        [str]
            The element symbols
        """
        return elements.to_symbols(self.asymmetric_atomic_numbers)

    @property
    def symbols(self):
        """The element symbols for the atoms in this configuration.

        Returns
        -------
        [str]
            The element symbols
        """
        if self.configuration.symmetry.n_symops == 1:
            return self.asymmetric_symbols
        else:
            return elements.to_symbols(self.atomic_numbers)

    @property
    def symmetry_matrices(self):
        """The 4x4 matrices for the symmetry operations."""
        return self.configuration.symmetry.symmetry_matrices

    @property
    def symops(self):
        """The symmetry operators as shorthand strings."""
        return self.configuration.symmetry.symops

    @symops.setter
    def symops(self, value):
        self.configuration.symmetry.symops = value

    @property
    def symop_to_atom(self):
        """List of list of symop #'s for creating symmetry atoms from asymmetric."""
        return self.configuration.symmetry.symop_to_atom

    @property
    def system_db(self):
        """The system database that we belong to."""
        return self._system_db

    @property
    def velocities(self):
        """The velocities as list of lists."""
        return self.get_velocities()

    @velocities.setter
    def velocities(self, xyz):
        """The velocities as list of lists."""
        return self.set_velocities(xyz)

    def add_attribute(
        self,
        name: str,
        coltype: str = "float",
        default: Any = None,
        notnull: bool = False,
        index: bool = False,
        pk: bool = False,
        references: str = None,
        on_delete: str = "cascade",
        on_update: str = "cascade",
        values: Any = None,
        configuration_dependent: bool = False,
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
                values=values,
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
                values=values,
            )

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more atoms.

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Need to handle the elements specially. Can give atomic numbers,
        # or symbols. By construction the references to elements are identical
        # to their atomic numbers.

        if "symbol" in kwargs:
            symbols = kwargs.pop("symbol")
            kwargs["atno"] = elements.to_atnos(symbols)

        # How many new rows there are
        n_rows, lengths = self._get_n_rows(**kwargs)

        # Fill in the atom table
        data = {}
        for column in self._atom_table.attributes:
            if column != "id" and column in kwargs:
                data[column] = kwargs.pop(column)

        if len(data) == 0:
            data["atno"] = [None] * n_rows

        ids = self._atom_table.append(n=n_rows, **data)

        # Now append to the coordinates table
        configuration = self.configuration.id
        data = {"configuration": configuration, "atom": ids}
        for column in self._coordinates_table.attributes:
            if column != "id" and column in kwargs:
                data[column] = kwargs.pop(column)
        self._coordinates_table.append(n=n_rows, **data)

        # And velocities table
        data = {"configuration": configuration, "atom": ids}
        have_velocities = False
        for column in self._velocities_table.attributes:
            if column != "id" and column in kwargs:
                data[column] = kwargs.pop(column)
                have_velocities = True
        if have_velocities:
            self._velocities_table.append(n=n_rows, **data)

        # And to the atomset
        table = _Table(self.system_db, "atomset_atom")
        table.append(atomset=self.atomset, atom=ids)

        self.configuration.symmetry.reset_atoms()

        return ids

    def atoms(self, *args):
        """Return an iterator over the atoms.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the atoms.
        """
        if self.have_velocities:
            columns = self._columns(velocities=True)
            column_defs = ", ".join(columns)

            # What tables are requested in the extra arguments?
            tables = set()
            if len(args) > 0:
                atom_columns = [*self._atom_table.attributes]
                coordinates_columns = [*self._coordinates_table.attributes]
                velocities_columns = [*self._velocities_table.attributes]
                for col, op, value in grouped(args, 3):
                    if "." in col:
                        tables.add(col.split(".")[0])
                    elif col in atom_columns:
                        tables.add("at")
                    elif col in coordinates_columns:
                        tables.add("co")
                    elif col in velocities_columns:
                        tables.add("ve")
                    else:
                        raise ValueError(f"Column '{col}' is not available")

            # Build the query based on the tables needed
            sql = (
                f"SELECT {column_defs}"
                " FROM atomset_atom AS aa, atom AS at, coordinates AS co, "
                "      velocities AS ve "
                "WHERE aa.atomset = ?"
                "  AND at.id = aa.atom"
                "  AND co.atom = at.id AND co.configuration = ?"
                "  AND ve.atom = at.id AND ve.configuration = ?"
            )
            parameters = [self.atomset, self.configuration.id, self.configuration.id]

            # And any extra selection criteria
            if len(args) > 0:
                for col, op, value in grouped(args, 3):
                    if op == "==":
                        op = "="
                    sql += f' AND "{col}" {op} ?'
                    parameters.append(value)
        else:
            columns = self._columns(velocities=False)
            column_defs = ", ".join(columns)

            # What tables are requested in the extra arguments?
            tables = set()
            if len(args) > 0:
                atom_columns = [*self._atom_table.attributes]
                coordinates_columns = [*self._coordinates_table.attributes]
                velocities_columns = [*self._velocities_table.attributes]
                for col, op, value in grouped(args, 3):
                    if "." in col:
                        tables.add(col.split(".")[0])
                    elif col in atom_columns:
                        tables.add("at")
                    elif col in coordinates_columns:
                        tables.add("co")
                    elif col in velocities_columns:
                        raise ValueError(
                            "Query for atom has velocities, but the atoms don't"
                        )
                    else:
                        raise ValueError(f"Column '{col}' is not available")

            # Build the query based on the tables needed
            sql = (
                f"SELECT {column_defs}"
                " FROM atomset_atom AS aa, atom AS at, coordinates AS co "
                "WHERE aa.atomset = ?"
                "  AND at.id = aa.atom"
                "  AND co.atom = at.id AND co.configuration = ?"
            )
            parameters = [self.atomset, self.configuration.id]

            # And any extra selection criteria
            if len(args) > 0:
                for col, op, value in grouped(args, 3):
                    if op == "==":
                        op = "="
                    sql += f' AND "{col}" {op} ?'
                    parameters.append(value)

        logger.debug("atoms query:")
        logger.debug(sql)
        logger.debug(parameters)
        logger.debug("---")

        return self.db.execute(sql, parameters)

    def diff(self, other):
        """Difference between these atoms and another

        Currently ignores velocities. Not sure what we want to do....

        Parameters
        ----------
        other : _Atoms
            The other atoms to diff against

        Result
        ------
        result : Dict
            The differences, described in a dictionary
        """
        result = {}

        # Check the columns
        columns = self._columns(velocities=False)
        other_columns = other._columns(velocities=False)

        column_defs = ", ".join(columns)
        other_column_defs = ", ".join(other_columns)

        if columns == other_columns:
            column_def = column_defs
        else:
            added = columns - other_columns
            if len(added) > 0:
                result["columns added"] = list(added)
            deleted = other_columns - columns
            if len(deleted) > 0:
                result["columns deleted"] = list(deleted)

            in_common = other_columns & columns
            if len(in_common) > 0:
                column_def = ", ".join(in_common)
            else:
                # No columns shared
                return result

        # Need to check the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        db = self.system_db
        other_db = other.system_db

        detach = False
        schema = self.schema
        if db.filename != other_db.filename:
            if db.is_attached(other_db):
                other_schema = db.attached_as(other_db)
            else:
                # Attach the other system_db in order to do comparisons.
                other_schema = self.system_db.attach(other_db)
                detach = True
        else:
            other_schema = other.schema

        atomset = self.atomset
        other_atomset = other.atomset

        changed = {}
        last = None
        sql = f"""
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {other_schema}.atom AS at, {other_schema}.coordinates AS co
           WHERE co.atom = at.id
             AND at.id
              IN (
                  SELECT atom
                    FROM {other_schema}.atomset_atom
                   WHERE atomset = {other_atomset}
                  )
          EXCEPT
          SELECT {column_def}
            FROM {schema}.atom AS at, {schema}.coordinates AS co
           WHERE co.atom = at.id
             AND at.id
              IN (
                  SELECT atom
                    FROM {schema}.atomset_atom
                    WHERE atomset = {atomset}
                 )
        )
         UNION ALL
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {schema}.atom AS at, {schema}.coordinates AS co
           WHERE co.atom = at.id
             AND at.id
              IN (
                  SELECT atom
                    FROM {schema}.atomset_atom
                   WHERE atomset = {atomset}
                 )
          EXCEPT
          SELECT {column_def}
            FROM {other_schema}.atom AS at, {other_schema}.coordinates AS co
           WHERE co.atom = at.id
             AND at.id
              IN (
                  SELECT atom
                    FROM {other_schema}.atomset_atom
                   WHERE atomset = {other_atomset}
                 )
        )
        ORDER BY id
        """

        for row in self.db.execute(sql):
            if last is None:
                last = row
            elif row["id"] == last["id"]:
                # changes = []
                changes = set()
                for k1, v1, v2 in zip(last.keys(), last, row):
                    if v1 != v2:
                        changes.add((k1, v1, v2))
                changed[row["id"]] = changes
                last = None
            else:
                last = row
        if len(changed) > 0:
            result["changed"] = changed

        # See about the rows added
        added = {}
        sql = f"""
        SELECT {column_defs}
          FROM {schema}.atom AS at, {schema}.coordinates AS co
         WHERE co.atom = at.id
           AND at.id
            IN (
                SELECT atom
                  FROM {schema}.atomset_atom
                 WHERE atomset = {atomset}
               )
           AND at.id
        NOT IN (
                SELECT atom
                  FROM {other_schema}.atomset_atom
                 WHERE atomset = {other_atomset}
               )
        """
        for row in self.db.execute(sql):
            added[row["id"]] = row[1:]

        if len(added) > 0:
            result["columns in added rows"] = row.keys()[1:]
            result["added"] = added

        # See about the rows deleted
        deleted = {}
        sql = f"""
        SELECT {other_column_defs}
          FROM {other_schema}.atom AS at, {other_schema}.coordinates AS co
         WHERE co.atom = at.id
           AND at.id
            IN (
                SELECT atom
                  FROM {other_schema}.atomset_atom
                 WHERE atomset = {other_atomset}
               )
           AND at.id
        NOT IN (
                SELECT atom
                  FROM {schema}.atomset_atom
                 WHERE atomset = {atomset}
               )
        """
        for row in self.db.execute(sql):
            deleted[row["id"]] = row[1:]

        if len(deleted) > 0:
            result["columns in deleted rows"] = row.keys()[1:]
            result["deleted"] = deleted

        # Detach the other database if needed
        if detach:
            self.system_db.detach(other_db)

        return result

    def get_as_dict(self, *args):
        """Return the atom data as a Python dictionary of lists.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        dict(str: [])
            A dictionary whose keys are the column names and values as lists
        """
        rows = self.atoms(*args)
        columns = [x[0] for x in rows.description]
        data = {key: [] for key in columns}
        for row in rows:
            for key, value in zip(columns, row):
                data[key].append(value)

        return data

    def get_ids(self, *args):
        """The ids of the selected atoms.

        Any extra arguments are triplets of column, operator, and value force
        additional selection criteria. The table names must be used in the column
        specification and are 'at' for the atom table and 'co' for the coordinate
        table.

        For example, if there are three arguments "at.atno", "=", "6" only the ids
        of the carbon atoms in the configuration will be returned.

        Parameters
        ----------
        args : [str]
            Further selection arguments, in sets of three:
            column, operator, and value, e.g. "at.atno = 6"

        Returns
        -------
        [int]
            The ids of the requested atoms.
        """

        # What tables are requested in the extra arguments?
        tables = set()
        if len(args) > 0:
            atom_columns = [*self._atom_table.attributes]
            coordinates_columns = [*self._coordinates_table.attributes]
            for col, op, value in grouped(args, 3):
                if "." in col:
                    tables.add(col.split(".")[0])
                elif col in atom_columns:
                    tables.add("at")
                elif col in coordinates_columns:
                    tables.add("co")
                else:
                    raise ValueError(f"Column '{col}' is not available")

        # Build the query based on the tables needed
        sql = "SELECT aa.atom FROM atomset_atom AS aa"
        if "at" in tables or "co" in tables:
            sql += ", atom AS at"
        if "co" in tables:
            sql += ", coordinates AS co"

        # The WHERE clauses bringing joining the tables
        sql += " WHERE aa.atomset = ?"
        parameters = [self.atomset]
        if "at" in tables or "co" in tables:
            sql += " AND at.id = aa.atom"
        if "co" in tables:
            sql += " AND co.atom = at.id AND co.configuration = ?"
            parameters.append(self.configuration.id)

        # And any extra selection criteria
        if len(args) > 0:
            for col, op, value in grouped(args, 3):
                if op == "==":
                    op = "="
                sql += f' AND "{col}" {op} ?'
                parameters.append(value)

        logger.debug("get_id query:")
        logger.debug(sql)
        logger.debug(parameters)
        logger.debug("---")

        return [x[0] for x in self.db.execute(sql, parameters)]

    def get_coordinates(
        self,
        fractionals=True,
        in_cell=False,
        as_array=False,
        asymmetric=False,
    ):
        """Return the coordinates optionally translated back into the principal
        unit cell.

        Parameters
        ----------
        fractionals : bool = True
            Return the coordinates as fractional coordinates for periodic
            systems. Non-periodic systems always use Cartesian coordinates.
        in_cell : bool, str = False
            Whether to translate the atoms into the unit cell, and if so
            whether to do so by molecule or just atoms.
        as_array : bool = False
            Whether to return the results as a np array or as a list of
            lists (the default).
        asymmetric : bool = False
            If true, return coordinates for only the symmetry-unique atoms.
            By default, expand to all atoms in the system.

        Returns
        -------
        abc : [N][float*3]
            The coordinates, either Cartesian or fractional
        """
        xyz = [[row["x"], row["y"], row["z"]] for row in self.atoms()]

        periodicity = self.configuration.periodicity
        if periodicity == 0:
            if asymmetric and self.configuration.symmetry.n_symops > 1:
                raise NotImplementedError("Point-group symmetry not handled yet.")
            if as_array:
                return np.array(xyz)
            else:
                return xyz

        cell = self.configuration.cell

        UVW = None
        if not asymmetric and self.n_symops > 1:
            # Get the asymmetric fractional coordinates as np array
            if self.configuration.coordinate_system == "Cartesian":
                UVW_asym = cell.to_fractionals(xyz, as_array=True)
            elif not isinstance(xyz, np.ndarray):
                UVW_asym = np.array(xyz)
            else:
                UVW_asym = xyz

            # Move into unit cell, remembering shift
            trans = np.floor(UVW_asym)
            UVW_asym -= trans

            # Apply the generating symmetry operators
            op = self.configuration.symmetry.symmetry_matrices
            generators = self.configuration.symmetry.atom_generators

            for i, indices in enumerate(generators):
                uvw4 = np.array([0.0, 0.0, 0.0, 1.0])
                uvw4[0:3] = UVW_asym[i]

                xformed = np.einsum("ijk,k", op[indices, :, :], uvw4)
                uvw = xformed[:, 0:3]
                if UVW is None:
                    UVW = uvw
                else:
                    UVW = np.concatenate((UVW, uvw))

        if (
            isinstance(in_cell, str)
            and "molecule" in in_cell
            and self.configuration.n_bonds > 0
        ):
            # Need fractionals...
            if UVW is None:
                if self.configuration.coordinate_system == "Cartesian":
                    UVW = cell.to_fractionals(xyz, as_array=True)
                elif not isinstance(xyz, np.ndarray):
                    UVW = np.array(xyz)
                else:
                    UVW = xyz

            molecules = self.configuration.find_molecules(as_indices=True)

            for indices in molecules:
                indices = np.array(indices)
                uvw_mol = np.take(UVW, indices, axis=0)
                center = np.average(uvw_mol, axis=0)
                delta = np.floor(center)
                uvw_mol -= delta
                np.put_along_axis(UVW, np.expand_dims(indices, axis=1), uvw_mol, axis=0)
            if fractionals:
                if as_array:
                    return UVW
                else:
                    return UVW.tolist()
            else:
                return cell.to_cartesians(UVW, as_array=as_array)
        elif in_cell:
            # Need fractionals...
            if UVW is None:
                if self.configuration.coordinate_system == "Cartesian":
                    UVW = cell.to_fractionals(xyz, as_array=True)
                elif not isinstance(xyz, np.ndarray):
                    UVW = np.array(xyz)
                else:
                    UVW = xyz
            delta = np.floor(UVW)
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
                if UVW is None:
                    if self.configuration.coordinate_system == "Cartesian":
                        return cell.to_fractionals(xyz, as_array=as_array)
                    elif as_array:
                        return np.array(xyz)
                    else:
                        return xyz
                else:
                    if as_array:
                        return UVW
                    else:
                        return UVW.tolist()
            else:
                if UVW is None:
                    if self.configuration.coordinate_system == "fractional":
                        return cell.to_cartesians(xyz, as_array=as_array)
                    elif as_array:
                        return np.array(xyz)
                    else:
                        return xyz
                else:
                    return cell.to_cartesians(UVW, as_array=as_array)

    def get_names(self, asymmetric=False):
        """Return the names of the atoms, return a default if not in the database.

        Parameters
        ----------
        asymmetric : bool = False
            Return just the names for the asymmetric atoms.

        Returns
        -------
        [str]
            The names of the atoms.
        """
        if "name" in self:
            name = self.get_column_data("name")
        else:
            name = []
            count = {}
            for symbol in self.asymmetric_symbols:
                if symbol not in count:
                    count[symbol] = 1
                else:
                    count[symbol] += 1
                name.append(f"{symbol}{count[symbol]}")

        symmetry = self.configuration.symmetry

        if asymmetric or symmetry.n_symops == 1:
            return name

        # Expand to the asymmetric atoms
        result = []
        count = {i: 0 for i in range(len(name))}
        for asym_atom in symmetry.atom_to_asymmetric_atom:
            count[asym_atom] += 1
            result.append(f"{name[asym_atom]}_{count[asym_atom]}")
        return result

    def get_n_atoms(self, *args):
        """Return the number of atoms meeting the cirteria.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        int
            The number of atoms matching the criteria.
        """
        # What tables are requested in the extra arguments?
        tables = set()
        if len(args) > 0:
            atom_columns = [*self._atom_table.attributes]
            coordinates_columns = [*self._coordinates_table.attributes]
            for col, op, value in grouped(args, 3):
                if "." in col:
                    tables.add(col.split(".")[0])
                elif col in atom_columns:
                    tables.add("at")
                elif col in coordinates_columns:
                    tables.add("co")
                else:
                    raise ValueError(f"Column '{col}' is not available")

        # Build the query based on the tables needed
        sql = "SELECT COUNT(*) FROM atomset_atom AS aa"
        if "at" in tables or "co" in tables:
            sql += ", atom AS at"
        if "co" in tables:
            sql += ", coordinates AS co"

        # The WHERE clauses bringing joining the tables
        sql += " WHERE aa.atomset = ?"
        parameters = [self.atomset]
        if "at" in tables or "co" in tables:
            sql += " AND at.id = aa.atom"
        if "co" in tables:
            sql += " AND co.atom = at.id AND co.configuration = ?"
            parameters.append(self.configuration.id)

        # And any extra selection criteria
        if len(args) > 0:
            for col, op, value in grouped(args, 3):
                if op == "==":
                    op = "="
                sql += f' AND "{col}" {op} ?'
                parameters.append(value)

        logger.debug("get_n_atoms query:")
        logger.debug(sql)
        logger.debug(parameters)
        logger.debug("---")

        self.cursor.execute(sql, parameters)
        return self.cursor.fetchone()[0]

    def get_velocities(
        self,
        fractionals=True,
        as_array=False,
    ):
        """Return the velocities.

        Symmetry is not supported, because it makes no (little?) sense for velocities.

        Parameters
        ----------
        fractionals : bool = True
            Return the velocities as fractional velocities for periodic
            systems. Non-periodic systems always use Cartesian velocities.
        as_array : bool = False
            Whether to return the results as a np array or as a list of
            lists (the default).

        Returns
        -------
        abc : [N][float*3]
            The velocities, either Cartesian or fractional
        """
        vxs = self.get_column_data("vx")
        vys = self.get_column_data("vy")
        vzs = self.get_column_data("vz")
        xyz = [[vx, vy, vz] for vx, vy, vz in zip(vxs, vys, vzs)]

        periodicity = self.configuration.periodicity
        if periodicity == 0:
            if as_array:
                return np.array(xyz)
            else:
                return xyz

        cell = self.configuration.cell

        if fractionals:
            if self.configuration.coordinate_system == "Cartesian":
                return cell.to_fractionals(xyz, as_array=as_array)
            elif as_array:
                return np.array(xyz)
            else:
                return xyz
        else:
            if self.configuration.coordinate_system == "fractional":
                return cell.to_cartesians(xyz, as_array=as_array)
            elif as_array:
                return np.array(xyz)
            else:
                return xyz

    def set_coordinates(self, xyz, fractionals=True):
        """Set the coordinates to new values.

        Parameters
        ----------
        fractionals : bool = True
            The coordinates are fractional coordinates for periodic
            systems. Ignored for non-periodic systems.

        Returns
        -------
        None
        """
        as_array = isinstance(xyz, np.ndarray)

        if as_array:
            n_coords = xyz.shape[0]
        else:
            n_coords = len(xyz)

        # May need to handle symmetry
        if n_coords != self.n_asymmetric_atoms and n_coords == self.n_atoms:
            xyz, error = self.configuration.symmetry.symmetrize_coordinates(
                xyz, fractionals=fractionals
            )

        x_column = self.get_column("x")
        y_column = self.get_column("y")
        z_column = self.get_column("z")

        xs = []
        ys = []
        zs = []

        periodicity = self.configuration.periodicity
        coordinate_system = self.configuration.coordinate_system
        if (
            periodicity == 0
            or (coordinate_system == "Cartesian" and not fractionals)
            or (coordinate_system == "fractional" and fractionals)
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
            cell = self.configuration.cell
            if coordinate_system == "fractional":
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

    def set_velocities(self, vxyz, fractionals=False):
        """Set the velocities to new values.

        Parameters
        ----------
        fractionals : bool = False
            The velocities are fractional velocities for periodic
            systems. Ignored for non-periodic systems.

        Returns
        -------
        None
        """
        if self.n_symops > 1:
            raise RuntimeError("Can't handle velocities with symmetry.")

        as_array = isinstance(vxyz, np.ndarray)

        vxs = []
        vys = []
        vzs = []

        periodicity = self.configuration.periodicity
        coordinate_system = self.configuration.coordinate_system
        if (
            periodicity == 0
            or (coordinate_system == "Cartesian" and not fractionals)
            or (coordinate_system == "fractional" and fractionals)
        ):
            if as_array:
                for vx, vy, vz in vxyz.tolist():
                    vxs.append(vx)
                    vys.append(vy)
                    vzs.append(vz)
            else:
                for vx, vy, vz in vxyz:
                    vxs.append(vx)
                    vys.append(vy)
                    vzs.append(vz)
        else:
            cell = self.configuration.cell
            if coordinate_system == "fractional":
                # Convert velocities to fractionals
                for vx, vy, vz in cell.to_fractionals(vxyz):
                    vxs.append(vx)
                    vys.append(vy)
                    vzs.append(vz)
            else:
                for vx, vy, vz in cell.to_cartesians(vxyz):
                    vxs.append(vx)
                    vys.append(vy)
                    vzs.append(vz)

        vx_column = self.get_column("vx")
        if len(vx_column) == 0:
            # No velocities in the database, so need to add rather than setFormatter
            self._velocities_table.append(
                n=len(vxs),
                vx=vxs,
                vy=vys,
                vz=vzs,
                atom=self.ids,
                configuration=self.configuration.id,
            )
        else:
            vy_column = self.get_column("vy")
            vz_column = self.get_column("vz")
            vx_column[0:] = vxs
            vy_column[0:] = vys
            vz_column[0:] = vzs

    def get_column(self, key: str) -> Any:
        """Get a Column object with the requested data

        Parameters
        ----------
        key : str
            The attribute to get.

        Returns
        -------
        Column
            A Column object containing the data.
        """
        if key in self._atom_table.attributes:
            sql = (
                f'SELECT at.rowid, at."{key}"'
                "   FROM atom as at, atomset_atom as aa"
                f" WHERE at.id = aa.atom AND aa.atomset = {self.atomset}"
            )
            return _Column(self._atom_table, key, sql=sql)
        elif key in self._coordinates_table.attributes:
            sql = (
                f'SELECT co.rowid, co."{key}"'
                "   FROM atom as at,"
                "        coordinates as co,"
                "        atomset_atom as aa"
                "  WHERE co.atom = at.id"
                f"   AND co.configuration = {self.configuration.id}"
                "    AND at.id = aa.atom"
                f"   AND aa.atomset = {self.atomset}"
            )
            return _Column(self._coordinates_table, key, sql=sql)
        elif key in self._velocities_table.attributes:
            sql = (
                f'SELECT ve.rowid, ve."{key}"'
                "   FROM atom as at,"
                "        velocities as ve,"
                "        atomset_atom as aa"
                "  WHERE ve.atom = at.id"
                f"   AND ve.configuration = {self.configuration.id}"
                "    AND at.id = aa.atom"
                f"   AND aa.atomset = {self.atomset}"
            )
            return _Column(self._velocities_table, key, sql=sql)
        else:
            raise KeyError(f"'{key}' not in atoms")

    def get_column_data(self, key: str) -> Any:
        """Return a column of data from the table.

        Parameters
        ----------
        key : str
            The attribute to get.

        Returns
        -------
        Column
            A Column object containing the data.
        """
        if key in self._atom_table.attributes:
            sql = (
                f'SELECT at."{key}"'
                "   FROM atom as at, atomset_atom as aa"
                f" WHERE at.id = aa.atom AND aa.atomset = {self.atomset}"
            )
            return [row[0] for row in self.db.execute(sql)]
        elif key in self._coordinates_table.attributes:
            sql = (
                f'SELECT co."{key}"'
                "   FROM atom as at,"
                "        coordinates as co,"
                "        atomset_atom as aa"
                "  WHERE co.atom = at.id"
                f"   AND co.configuration = {self.configuration.id}"
                "    AND at.id = aa.atom"
                f"   AND aa.atomset = {self.atomset}"
            )
            return [row[0] for row in self.db.execute(sql)]
        elif key in self._velocities_table.attributes:
            sql = (
                f'SELECT ve."{key}"'
                "   FROM atom as at,"
                "        velocities as ve,"
                "        atomset_atom as aa"
                "  WHERE ve.atom = at.id"
                f"   AND ve.configuration = {self.configuration.id}"
                "    AND at.id = aa.atom"
                f"   AND aa.atomset = {self.atomset}"
            )
            return [row[0] for row in self.db.execute(sql)]
        else:
            raise KeyError(f"'{key}' not in atoms")

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
                        'key "{}" has the wrong number of values, '.format(key)
                        + "{}. Should be 1 or the number of atoms ({}).".format(
                            length, n_rows
                        )
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

    def delete(self, atoms) -> int:
        """Delete the atoms listed

        Parameters
        ----------
        atoms : [int]
            The list of atoms to delete, or 'all' or '*'

        Returns
        -------
        None
        """
        # Delete the listed atoms, which will cascade to delete coordinates and
        # velocities
        if atoms == "all" or atoms == "*":
            sql = """
            DELETE FROM atom
             WHERE id IN (SELECT atom FROM atomset_atom WHERE atomset = ?)
            """
            parameters = (self.atomset,)
            self.db.execute(sql, parameters)
        else:
            sql = """
            DELETE FROM atom
             WHERE id = ?
               AND id IN (SELECT atom FROM atomset_atom WHERE atomset = ?)
            """
            parameters = [(i, self.atomset) for i in atoms]
            self.db.executemany(sql, parameters)

        self.configuration.symmetry.reset_atoms()

    def to_dataframe(self):
        """Return the contents of the table as a Pandas Dataframe."""
        data = {}
        rows = self.atoms()
        for row in rows:
            data[row[0]] = row[1:]

        columns = [x[0] for x in rows.description[1:]]

        df = pandas.DataFrame.from_dict(data, orient="index", columns=columns)

        return df

    def _columns(self, velocities=True):
        """The list of columns across the atom, coordinates and velocities tables.

        Uses 'at', 'co', and 've' as the shorthand for the full table names.
        """
        atom_columns = [*self._atom_table.attributes]
        coordinates_columns = [*self._coordinates_table.attributes]
        coordinates_columns.remove("atom")
        columns = [f'at."{x}"' for x in atom_columns]
        columns += [f'co."{x}"' for x in coordinates_columns]
        if velocities:
            velocities_columns = [*self._velocities_table.attributes]
            velocities_columns.remove("atom")
            columns += [f've."{x}"' for x in velocities_columns]

        return columns


class _SubsetAtoms(_Atoms):
    """The atoms in a subset.

    Parameters
    ----------
    configuration : _Configuration
        The configuration of interest.
    subset_id : int
        The id of the subset.
    template_order : bool
        Whether to return atoms and properties in the order of the template
        if the template is full. Defaults to True.
    """

    def __init__(self, configuration, subset_id) -> None:
        self._sid = subset_id
        self.template_order = True

        # Caching
        self._template_id = None
        self._template = None

        super().__init__(configuration)

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError("Not implemented for subsets, yet!")

    @property
    def subset_id(self):
        """The subset for these atoms."""
        return self._sid

    @property
    def atomic_numbers(self):
        """The atomic numbers of the subset atoms.

        Note that subsets refer to asymmetric atoms!

        Returns
        -------
        [int]
            The atomic numbers.
        """
        return self.get_column_data("atno")

    @property
    def atomic_masses(self):
        """The atomic masses of the atoms in the subset.

        Note that subsets refer to asymmetric atoms!

        Returns
        -------
        [int]
            The atomic numbers.
        """

        if "mass" in self:
            result = self.get_column_data("mass")
        else:
            atnos = self.atomic_numbers
            result = elements.masses(atnos)
        return result

    @property
    def have_velocities(self):
        """Whether there are any velocities for this configuration."""
        sql = (
            "SELECT COUNT(*)"
            " FROM subset_atom AS sa, velocities AS ve "
            "WHERE sa.subset = ?"
            "  AND ve.atom = sa.atom AND ve.configuration = ?"
        )
        parameters = [self.subset_id, self.configuration.id]
        self.cursor.execute(sql, parameters)
        return self.cursor.fetchone()[0] > 0

    @property
    def ids(self):
        """The ids of atoms in this subset."""
        sql = """
        SELECT atom
          FROM subset_atom
         WHERE subset = ?
        """
        if self.template.is_full and self.template_order:
            sql += "ORDER BY templateatom"

        return [x[0] for x in self.db.execute(sql, (self.subset_id,))]

    @property
    def n_atoms(self) -> int:
        """The number of atoms in this subset."""
        sql = """
        SELECT COUNT(*)
          FROM atomset_atom
         WHERE atomset = ?
           AND atom IN (SELECT atom FROM subset_atom WHERE subset = ?)
        """
        self.cursor.execute(sql, (self.atomset, self.subset_id))
        return self.cursor.fetchone()[0]

    @property
    def template(self):
        """The template for this subset."""
        if self._template is None:
            self._template = self.system_db.templates.get(self.template_id)
        return self._template

    @property
    def template_id(self):
        """The id of the template for this subset."""
        if self._template_id is None:
            sql = "SELECT template FROM SUBSET WHERE id = ?"
            self.cursor.execute(sql, (self.subset_id,))
            self._template_id = self.cursor.fetchone()[0]
        return self._template_id

    def add(self, ids):
        """Add atoms to the subset.

        Parameters
        ----------
        ids : [int]
            The ids of the atoms to add. They are silently ignored if
            they are already in the subset.

        Raises
        ------
        ValueError
            If this subset has a full template or if an atom is not in the
            configuration.
        """
        # Check if the template has a template configuration
        if self.template.is_full:
            raise ValueError("Cannot add atoms to a subset for a full template.")

        # Remove any ids already in the subset
        atom_ids = set(self.ids)
        for aid in atom_ids & set(ids):
            ids.remove(aid)

        sa = self.system_db["subset_atom"]
        sa.append(subset=self.subset_id, atom=ids)

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more atoms.

        This is not allowed for a subset.

        Raises
        ------
        RuntimeError
            Adding atoms to a configuration is not allowed from a subset.
        """
        raise RuntimeError(
            "Adding atoms to a configuration is not allowed from a subset."
        )

    def atoms(self, *args):
        """Return an iterator over the atoms.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the atoms.
        """
        if self.have_velocities:
            columns = self._columns()
            column_defs = ", ".join(columns)

            sql = f"""
            SELECT {column_defs}
              FROM atom as at, coordinates as co, velocities as ve, subset_atom as sa
             WHERE co.atom = at.id
               AND co.configuration = ?
               AND ve.atom = at.id
               AND ve.configuration = ?
               AND at.id = sa.atom
               AND sa.subset = ?
            """

            parameters = [self.configuration.id, self.configuration.id, self.subset_id]
            if len(args) > 0:
                for col, op, value in grouped(args, 3):
                    if op == "==":
                        op = "="
                    sql += f' AND "{col}" {op} ?'
                    parameters.append(value)
        else:
            columns = self._columns(velocities=False)
            column_defs = ", ".join(columns)

            sql = f"""
            SELECT {column_defs}
              FROM atom as at, coordinates as co, subset_atom as sa
             WHERE co.atom = at.id
               AND co.configuration = ?
               AND at.id = sa.atom
               AND sa.subset = ?
            """

            parameters = [self.configuration.id, self.subset_id]
            if len(args) > 0:
                for col, op, value in grouped(args, 3):
                    if op == "==":
                        op = "="
                    sql += f' AND "{col}" {op} ?'
                    parameters.append(value)

            if self.template.is_full and self.template_order:
                sql += "ORDER BY sa.templateatom"

        return self.db.execute(sql, parameters)

    def get_n_atoms(self, *args):
        """Return the number of atoms meeting the criteria.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        int
            The number of atoms matching the criteria.
        """
        sql = """
        SELECT COUNT(*)
          FROM atom as at, coordinates as co, subset_atom as sa
         WHERE co.atom = at.id
           AND co.configuration = ?
           AND at.id = sa.atom
           AND sa.subset = ?
        """
        parameters = [self.configuration.id, self.subset_id]
        if len(args) > 0:
            for col, op, value in grouped(args, 3):
                if op == "==":
                    op = "="
                sql += f' AND "{col}" {op} ?'
                parameters.append(value)

        self.cursor.execute(sql, parameters)
        return self.cursor.fetchone()[0]

    def diff(self, other):
        """Difference between these atoms and another

        Parameters
        ----------
        other : _Atoms
            The other atoms to diff against

        Result
        ------
        result : Dict
            The differences, described in a dictionary
        """
        raise NotImplementedError()

    def get_column(self, key: str) -> Any:
        """Get a Column object with the requested data

        Parameters
        ----------
        key : str
            The attribute to get.

        Returns
        -------
        Column
            A Column object containing the data.
        """
        if key in self._atom_table.attributes:
            sql = f"""
            SELECT at.rowid, at."{key}"
              FROM atom as at, subset_atom as sa
             WHERE at.id = sa.atom
               AND sa.subset = {self.subset_id}
            """
            if self.template.is_full and self.template_order:
                sql += "ORDER BY sa.templateatom"
            return _Column(self._atom_table, key, sql=sql)
        elif key in self._coordinates_table.attributes:
            sql = f"""
            SELECT co.rowid, co."{key}"
              FROM coordinates as co, subset_atom as sa
             WHERE co.atom = sa.atom
               AND co.configuration = {self.configuration.id}
               AND sa.subset = {self.subset_id}
            """
            if self.template.is_full and self.template_order:
                sql += "ORDER BY sa.templateatom"
            return _Column(self._coordinates_table, key, sql=sql)
        else:
            raise KeyError(f"'{key}' not in atoms")

    def get_column_data(self, key: str) -> Any:
        """Return a column of data from the table.

        Parameters
        ----------
        key : str
            The attribute to get.

        Returns
        -------
        Column
            A Column object containing the data.
        """
        if key in self._atom_table.attributes:
            sql = f"""
            SELECT at."{key}"
              FROM atom as at, subset_atom as sa
             WHERE at.id = sa.atom
               AND sa.subset = {self.subset_id}
            """
            if self.template.is_full and self.template_order:
                sql += "ORDER BY sa.templateatom"
            return [row[0] for row in self.db.execute(sql)]
        elif key in self._coordinates_table.attributes:
            sql = f"""
            SELECT co."{key}"
              FROM coordinates as co, subset_atom as sa
             WHERE co.atom = sa.atom
               AND co.configuration = {self.configuration.id}
               AND sa.subset = {self.subset_id}
            """
            if self.template.is_full and self.template_order:
                sql += "ORDER BY sa.templateatom"
            return [row[0] for row in self.db.execute(sql)]
        else:
            raise KeyError(f"'{key}' not in atoms")

    def get_ids(self, *args):
        """The ids of the atoms.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        [int]
            The ids of the requested atoms.
        """

        sql = """
        SELECT at.id
          FROM atom as at, coordinates as co, subset_atom as sa
         WHERE at.id = sa.atom
           AND co.configuration = ?
           AND sa.subset = ?
        """

        parameters = [self.configuration.id, self.subset_id]
        if len(args) > 0:
            for col, op, value in grouped(args, 3):
                if op == "==":
                    op = "="
                sql += f' AND "{col}" {op} ?'
                parameters.append(value)

        if self.template.is_full and self.template_order:
            sql += "ORDER BY sa.templateatom"

        return [x[0] for x in self.db.execute(sql, parameters)]

    def delete(self, atoms) -> int:
        """Delete the atoms listed

        Parameters
        ----------
        atoms : [int]
            The list of atoms to delete, or 'all' or '*'

        Returns
        -------
        None
        """
        raise NotImplementedError()

    def remove(self, ids):
        """Remove the given atoms from the subset.

        Parameters
        ----------
        ids : [int]
            The ids of the atoms to delete.

        Raises
        ------
        ValueError
            If this subset has a full template or if an atom is not in the
            subset.
        """
        # Check if the template has a template configuration
        if self.template.is_full:
            raise ValueError("Cannot add atoms to a subset for a full template.")

        # Check that the atoms are in this subset!
        atom_ids = self.ids
        for aid in ids:
            if aid not in atom_ids:
                raise ValueError(f"Atom id={aid} is not in the subset.")

        parameters = [(_id, self.subset_id) for _id in ids]
        self.db.executemany(
            "DELETE FROM subset_atom WHERE atom=? AND subset=?", parameters
        )

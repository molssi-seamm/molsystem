# -*- coding: utf-8 -*-

"""A dictionary-like object for holding bonds

Based on tables in an SQLite database.
"""

from collections.abc import Sequence
from itertools import zip_longest
import logging
import sqlite3

import numpy as np
import pandas

from .column import _Column
from .table import _Table
from .frozencolumn import _FrozenColumn

logger = logging.getLogger(__name__)


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Bonds(_Table):
    """The Bonds class describes the bonds in the system.

    :meta public:
    """

    def __init__(self, configuration):

        self._configuration = configuration

        self._system = None

        self._bondset = self._configuration.bondset

        super().__init__(configuration._system_db, "bond")

    def __enter__(self):
        """Copy the tables to a backup for a 'with' statement."""
        self.system_db["bondset_bond"].__enter__()
        self.system_db["bond"].__enter__()
        return self

    def __exit__(self, etype, value, traceback):
        """Handle returning from a 'with' statement."""
        if etype is None:
            self.configuration.version = self.configuration.version + 1
        self.system_db["bondset_bond"].__exit__(etype, value, traceback)
        return self.system_db["bond"].__exit__(etype, value, traceback)

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    @property
    def ids(self):
        """The ids of the bonds."""
        return self.get_ids()

    @property
    def bondorders(self):
        """The bond orders."""
        return self.get_column_data("bondorder")

    @property
    def bondset(self):
        """The bondset for these bonds."""
        return self._bondset

    @property
    def configuration(self):
        """Return the configuration."""
        return self._configuration

    @property
    def cursor(self):
        """The database connection."""
        return self.system_db.cursor

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def n_bonds(self):
        """The number of bonds."""
        self.cursor.execute(
            "SELECT COUNT(*) FROM bondset_bond WHERE bondset = ?", (self.bondset,)
        )
        return self.cursor.fetchone()[0]

    @property
    def system_db(self):
        """The system database that we belong to."""
        return self._system_db

    @property
    def system(self):
        """The system that we belong to."""
        if self._system is None:
            self.cursor.execute(
                "SELECT system FROM configuration WHERE id = ?", (self.id,)
            )
            self._system = self.system_db.get_system(self.cursor.fetchone()[0])
        return self._system

    def append(self, **kwargs):
        """Append one or more bonds

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.

        Parameters
        ----------
        bonds : [sqlite3.Rows] or {str : [int]}
            The bonds as either SQLite Rows or a dictionary of values or lists.
            Any other arguments override the same value in 'bonds'.
        i : [int]
            The atom indices of the first atom.
        j : [int]
            The atom indices of the second atom.
        bondorder : [int] optional
            Bond orders, defaults to 1 for single bonds.
        other : [int, float, str] optional
            Other optional attributes of bonds.
        """
        if "bonds" in kwargs:
            bonds = kwargs.pop("bonds")
            if isinstance(bonds, sqlite3.Row):
                for item in bonds.keys():
                    if item != "id":
                        kwargs[item] = bonds[item]
            elif isinstance(bonds, Sequence) and isinstance(bonds[0], sqlite3.Row):
                for item in bonds[0].keys():
                    if item != "id":
                        kwargs[item] = [row[item] for row in bonds]
            else:
                try:
                    for key, value in bonds.items():
                        if isinstance(value, Sequence):
                            kwargs[key] = value
                        else:
                            kwargs[key] = [value]
                except AttributeError:
                    raise TypeError(
                        "bonds argument must be sqlite3.Row or a mapping "
                        f"(dict) type, not {type(bonds)}."
                    )

        # Check keys and lengths of added bonds
        if "i" not in kwargs or "j" not in kwargs:
            raise KeyError("The atoms i & j are required!")

        i = kwargs.pop("i")
        j = kwargs.pop("j")
        try:
            len_i = len(i)
        except TypeError:
            len_i = 1
            i = (i,)
        try:
            len_j = len(j)
        except TypeError:
            len_j = 1
            j = (j,)

        if len_i == 0 and len_j == 0:
            return

        if len_i == 1:
            if len_j > 1:
                i = len_j * [i[0]]
        elif len_j == 1:
            j = len_i * [j[0]]
        elif len_i != len_j:
            raise IndexError(
                f'key "j" has the wrong number of values, {len_j}. '
                f"Should be 1 or the number of values in i, {len_i}."
            )
        # Ensure that i < j

        if isinstance(i, np.ndarray):
            i = i.tolist()
        if isinstance(j, np.ndarray):
            j = j.tolist()

        i2 = []
        j2 = []
        for i_, j_ in zip(i, j):
            # will need to handle offsets here at some point
            if not isinstance(i_, int) or not isinstance(j_, int):
                raise TypeError("'i' and 'j', the atom indices, must be integers")
            if i_ < j_:
                i2.append(i_)
                j2.append(j_)
            else:
                i2.append(j_)
                j2.append(i_)

        ids = super().append(i=i2, j=j2, **kwargs)

        # And to the bondset
        table = _Table(self.system_db, "bondset_bond")
        table.append(bondset=self.bondset, bond=ids)

        return ids

    def bonds(self, *args):
        """Returns an iterator over the bonds.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the bonds.
        """
        sql = (
            "SELECT bond.* FROM bond, bondset_bond"
            " WHERE bond.id = bondset_bond.bond"
            "   AND bondset_bond.bondset = ?"
        )
        return self.db.execute(sql, (self.bondset,))

    def contains_bond(self, i, j):
        """Whether there is a bond between atoms i and j."""
        # get canonical order
        if i > j:
            j, i = i, j
        sql = (
            "SELECT COUNT(*) FROM bond, bondset_bond"
            " WHERE bond.i = ? AND bond.j = ?"
            "   AND bondset_bond.bondset = ?"
        )
        self.cursor.execute(sql, (i, j, self.bondset))
        return self.cursor.fetchone()[0] == 1

    def delete_bond(self, i, j=None, force=False):
        """Delete a bond, if present.

        Parameters
        ----------
        i : int
            Either the first atom in the bond or the bond id.
        j : int = None
            If not None, then the second atom in the bond.
        force : bool = False
            If true, ignore missing bonds.
        """
        # get canonical order
        if j is None:
            bond_id = i
        else:
            if i > j:
                j, i = i, j
            sql = (
                "SELECT bond.id FROM bond, bondset_bond"
                " WHERE bond.i = ? AND bond.j = ?"
                "   AND bondset_bond.bondset = ?"
            )
            self.cursor.execute(sql, (i, j, self.bondset))
            bond_id = self.cursor.fetchone()
            if bond_id is not None:
                bond_id = bond_id[0]

        if bond_id is None:
            if not force:
                raise ValueError("The bond did not exist.")
        else:
            self.cursor.execute("DELETE FROM bond WHERE id = ?", (bond_id,))

    def delete(self, atoms=None):
        """Deletes all the bonds, optionally from atoms.

        Parameters
        ----------
        atoms : [int] = None
            The list of atoms whose bonds are deleted

        Returns
        -------
        None
        """
        if atoms is not None:
            sql = (
                "DELETE FROM bond"
                " WHERE (i = ? OR j = ?)"
                "   AND id in ("
                "              SELECT bond FROM bondset_bond"
                "              WHERE bondset = ?"
                "             )"
            )
            parameters = [(i, i, self.bondset) for i in atoms]
            self.db.executemany(sql, parameters)
        else:
            sql = (
                "DELETE FROM bond"
                " WHERE id in ("
                "              SELECT bond FROM bondset_bond"
                "               WHERE bondset = ?"
                "             )"
            )
            self.db.execute(sql, (self.bondset,))

    def diff(self, other):
        """Difference between these bonds and another

        Parameters
        ----------
        other : _Bonds
            The other bonds to diff against

        Result
        ------
        result : Dict
            The differences, described in a dictionary
        """
        result = {}

        # Check the columns
        columns = self._columns()
        other_columns = other._columns()

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

        bondset = self.bondset
        other_bondset = other.bondset

        changed = {}
        last = None
        sql = f"""
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {other_schema}.bond
           WHERE id
              IN (
                  SELECT bond
                    FROM {other_schema}.bondset_bond
                   WHERE bondset = {other_bondset}
                  )
          EXCEPT
          SELECT {column_def}
            FROM {schema}.bond
           WHERE id
              IN (
                  SELECT bond
                    FROM {schema}.bondset_bond
                    WHERE bondset = {bondset}
                 )
        )
         UNION ALL
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {schema}.bond
           WHERE id
              IN (
                  SELECT bond
                    FROM {schema}.bondset_bond
                   WHERE bondset = {bondset}
                 )
          EXCEPT
          SELECT {column_def}
            FROM {other_schema}.bond
           WHERE id
              IN (
                  SELECT bond
                    FROM {other_schema}.bondset_bond
                   WHERE bondset = {other_bondset}
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
          FROM {schema}.bond
         WHERE id
            IN (
                SELECT bond
                  FROM {schema}.bondset_bond
                 WHERE bondset = {bondset}
               )
           AND id
        NOT IN (
                SELECT bond
                  FROM {other_schema}.bondset_bond
                 WHERE bondset = {other_bondset}
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
          FROM {other_schema}.bond
         WHERE id
            IN (
                SELECT bond
                  FROM {other_schema}.bondset_bond
                 WHERE bondset = {other_bondset}
               )
           AND id
        NOT IN (
                SELECT bond
                  FROM {schema}.bondset_bond
                 WHERE bondset = {bondset}
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
        """Return the bond data as a Python dictionary of lists.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        dict(str: [])
            A dictionary whose keys are the column names and values as lists
        """
        rows = self.bonds(*args)
        columns = [x[0] for x in rows.description]
        data = {key: [] for key in columns}
        for row in rows:
            for key, value in zip(columns, row):
                data[key].append(value)

        return data

    def get_bond(self, i, j=None, force=False):
        """Return the row for a bond.

        Parameters
        ----------
        i : int
            Either the first atom in the bond or the bond id.
        j : int = None
            If not None, then the second atom in the bond.
        force : bool = False
            If true, ignore missing bonds.

        Returns
        -------
        row : SQLite3.Row
            The Row object for the bond
        """
        if j is None:
            sql = (
                "SELECT bond.* FROM bond, bondset_bond"
                " WHERE id = ?"
                "   AND bondset_bond.bond = bond.id"
                "   AND bondset_bond.bondset = ?"
            )
            self.cursor.execute(sql, (i, self.bondset))
            row = self.cursor.fetchone()
            if row is None and not force:
                raise KeyError(f"No bond id = {i} found")
        else:
            if i > j:
                j, i = i, j
            sql = (
                "SELECT bond.* FROM bond, bondset_bond"
                " WHERE bond.i = ? AND bond.j = ?"
                "   AND bondset_bond.bondset = ?"
            )
            self.cursor.execute(sql, (i, j, self.bondset))
            row = self.cursor.fetchone()
            if row is None and not force:
                raise KeyError(f"No bond from {i} to {j} found")
        return row

    def get_column(self, key):
        """Return a column of data.

        Parameters
        ----------
        key : str
            The column (attribute) to getLogger

        Returns
        -------
        column : Column
            The column object requested.
        """
        if key == "i" or key == "j" or key == "id":
            sql = f"""
            SELECT {key}
               FROM bond, bondset_bond
              WHERE bond.id = bondset_bond.bond
                AND bondset_bond.bondset = {self.bondset}
            """
            return _FrozenColumn(self, key, sql)
        else:
            sql = f"""
            SELECT id, {key}
               FROM bond, bondset_bond
              WHERE bond.id = bondset_bond.bond
                AND bondset_bond.bondset = {self.bondset}
            """
            return _Column(self, key, sql)

    def get_column_data(self, key):
        """Return a column of data as a Python list.

        Parameters
        ----------
        key : str
            The column (attribute) to getLogger

        Returns
        -------
        column : Column
            The column object requested.
        """
        sql = f"""
        SELECT {key}
          FROM bond, bondset_bond
         WHERE bond.id = bondset_bond.bond
           AND bondset_bond.bondset = {self.bondset}
        """
        return [row[0] for row in self.db.execute(sql)]

    def get_ids(self, *args):
        """The ids of the bonds.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        [int]
            The ids of the requested bonds.
        """

        sql = (
            "SELECT bond.id FROM bond, bondset_bond"
            " WHERE bond.id = bondset_bond.bond"
            "   AND bondset_bond.bondset = ?"
        )

        parameters = [self.bondset]
        if len(args) > 0:
            for col, op, value in grouped(args, 3):
                if op == "==":
                    op = "="
                sql += f' AND "{col}" {op} ?'
                parameters.append(value)

        return [x[0] for x in self.db.execute(sql, parameters)]

    def import_bonds(self, ids):
        """Import existing bonds into this configuration.

        Parameters
        ----------
        ids : iterable(int)
            The ids of the bonds to import.
        """
        table = _Table(self.system_db, "bondset_bond")
        table.append(bondset=self.bondset, bond=ids)

    def to_dataframe(self):
        """Return the bonds as a Pandas Dataframe."""
        cursor = self.bonds()
        data = {row[0]: row[1:] for row in cursor}
        columns = [x[0] for x in cursor.description[1:]]

        df = pandas.DataFrame.from_dict(data, orient="index", columns=columns)

        return df


class _SubsetBonds(_Bonds):
    """The Bonds class describes the bonds in a subset

    :meta public:
    """

    def __init__(self, configuration, subset_id):
        self._sid = subset_id

        # Caching
        self._template_id = None
        self._template = None

        super().__init__(configuration)

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    @property
    def ids(self):
        """The ids of the bonds in the subset."""
        return self.get_bond_ids()

    @property
    def n_bonds(self):
        """The number of bonds."""
        sql = """
        SELECT COUNT(*) FROM bond
         WHERE id IN (SELECT bond FROM bondset_bond WHERE bondset = ?)
           AND i IN (SELECT atom FROM subset_atom WHERE subset = ?)
           AND j IN (SELECT atom FROM subset_atom WHERE subset = ?)
        """
        self.cursor.execute(sql, (self.bondset, self.subset_id))
        return self.cursor.fetchone()[0]

    def append(self, **kwargs):
        """Append one or more bonds.

        Not currently allowed in subsets."""
        raise NotImplementedError("Can't add bonds in a subset yet.")

    def bonds(self, *args):
        """Returns an iterator over the bonds.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the bonds.
        """
        sql = """
        SELECT * FROM bond
         WHERE id IN (SELECT bond FROM bondset_bond WHERE bondset = ?)
           AND i IN (SELECT atom FROM subset_atom WHERE subset = ?)
           AND j IN (SELECT atom FROM subset_atom WHERE subset = ?)
        """
        return self.db.execute(sql, (self.bondset, self.subset_id))

    def contains_bond(self, i, j):
        """Whether there is a bond between atoms i and j."""
        # get canonical order
        if i > j:
            j, i = i, j
        sql = """
        SELECT COUNT(*) FROM bond
         WHERE i = ? AND j = ?
           AND id IN (SELECT bond FROM bondset_bond WHERE bondset = ?)
           AND i IN (SELECT atom FROM subset_atom WHERE subset = ?)
           AND j IN (SELECT atom FROM subset_atom WHERE subset = ?)
        """
        self.cursor.execute(sql, (i, j, self.bondset, self.subset_id))
        return self.cursor.fetchone()[0] == 1

    def delete_bond(self, i, j=None, force=False):
        """Delete a bond, if present.


        Not currently allowed in subsets."""
        raise NotImplementedError("Can't delete bonds from a subset yet.")

    def delete(self, atoms=None):
        """Deletes all the bonds, optionally from atoms.

        Not currently allowed in subsets."""
        raise NotImplementedError("Can't add bonds in a subset yet.")

    def diff(self, other):
        """Difference between these bonds and another

        Not currently allowed in subsets."""
        raise NotImplementedError("Can't add bonds in a subset yet.")

    def get_bond(self, i, j=None, force=False):
        """Return the row for a bond.

        Parameters
        ----------
        i : int
            Either the first atom in the bond or the bond id.
        j : int = None
            If not None, then the second atom in the bond.
        force : bool = False
            If true, ignore missing bonds.

        Returns
        -------
        row : SQLite3.Row
            The Row object for the bond
        """
        if j is None:
            sql = """
            SELECT * FROM bond
             WHERE id = ?
               AND id IN (SELECT bond FROM bondset_bond WHERE bondset = ?)
               AND i IN (SELECT atom FROM subset_atom WHERE subset = ?)
               AND j IN (SELECT atom FROM subset_atom WHERE subset = ?)
            """
            self.cursor.execute(sql, (i, self.bondset, self.subset_id, self.subset_id))
            row = self.cursor.fetchone()
            if row is None and not force:
                raise KeyError(f"No bond id = {i} found")
        else:
            if i > j:
                j, i = i, j
            sql = """
            SELECT COUNT(*) FROM bond
             WHERE i = ? AND j = ?
               AND id IN (SELECT bond FROM bondset_bond WHERE bondset = ?)
               AND i IN (SELECT atom FROM subset_atom WHERE subset = ?)
               AND j IN (SELECT atom FROM subset_atom WHERE subset = ?)
            """
            self.cursor.execute(
                sql, (i, j, self.bondset, self.subset_id, self.subset_id)
            )
            row = self.cursor.fetchone()
            if row is None and not force:
                raise KeyError(f"No bond from {i} to {j} found")
        return row

    def get_column(self, key):
        """Return a column of data.

        Parameters
        ----------
        key : str
            The column (attribute) to getLogger

        Returns
        -------
        column : Column
            The column object requested.
        """
        if key == "i" or key == "j" or key == "id":
            sql = f"""
            SELECT {key} FROM bond
               AND id IN (
                  SELECT bond FROM bondset_bond WHERE bondset = {self.bondset}
               )
               AND i IN (
                  SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
               )
               AND j IN (
                  SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
            )
            """
            return _FrozenColumn(self, key, sql)
        else:
            sql = f"""
            SELECT id, {key} FROM bond
               AND id IN (
                  SELECT bond FROM bondset_bond WHERE bondset = {self.bondset}
               )
               AND i IN (
                  SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
               )
               AND j IN (
                  SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
            )
            """
            return _Column(self, key, sql)

    def get_column_data(self, key):
        """Return a column of data as a Python list.

        Parameters
        ----------
        key : str
            The column (attribute) to getLogger

        Returns
        -------
        column : Column
            The column object requested.
        """
        sql = f"""
        SELECT {key} FROM bond
           AND id IN (
              SELECT bond FROM bondset_bond WHERE bondset = {self.bondset}
           )
           AND i IN (
              SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
           )
           AND j IN (
              SELECT atom FROM subset_atom WHERE subset = {self.subset_id}
        )
        """
        return [row[0] for row in self.db.execute(sql)]

# -*- coding: utf-8 -*-

"""A dictionary-like object for holding a system
"""

import collections.abc
import logging
import os
import shutil
import sqlite3
import tempfile
from typing import Any, Dict

# import molsystem
from molsystem.table import _Table
from molsystem.atoms import _Atoms
import seamm_util

logger = logging.getLogger(__name__)


class System(collections.abc.MutableMapping):

    def __init__(self, other=None, **kwargs):
        # these are the only things not copied
        self._checkpoints = []
        self._filename = None
        self._db = None
        self._cursor = None

        # copy constructor
        if other and isinstance(other, System):
            raise NotImplementedError(
                'copy constructor is under construction!'
            )
        else:
            pass

        if 'filename' in kwargs:
            self.filename = kwargs.pop('filename')
        else:
            self.filename = 'system.seamm'

    def __enter__(self) -> Any:
        self.db.commit()

        tmpdir = tempfile.mkdtemp()
        filepath = os.path.join(tmpdir, 'backup.seamm')
        backup = sqlite3.connect(filepath)
        backup.row_factory = sqlite3.Row
        self._db.execute('PRAGMA foreign_keys = ON')

        self.db.backup(backup)

        self._checkpoints.append((backup, tmpdir, filepath))
        return self

    def __exit__(self, etype, value, traceback) -> None:
        backup, tmpdir, backup_file = self._checkpoints.pop()

        if etype is None:
            self.db.commit()

            # Log the changes
            self._log_changes(backup, backup_file)

            # Not sure why this commit is needed...
            self.db.commit()
        else:
            self.cursor.close()
            self.db.close()
            self._db = sqlite3.connect(self._filename)
            self._db.row_factory = sqlite3.Row
            self._db.execute('PRAGMA foreign_keys = ON')
            self._cursor = self._db.cursor()
            backup.backup(self._db)

        # and delete the copy
        backup.close()
        shutil.rmtree(tmpdir)

    def __getitem__(self, key):
        """Allow [] access to the dictionary!"""
        if key == 'atoms':
            return _Atoms(self, key)
        else:
            return _Table(self, key)

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        raise NotImplementedError(f"Table '{key}' cannot be created yet")

    def __delitem__(self, key):
        """Allow deletion of keys"""
        if key in self:
            self.cursor.execute(f"DROP TABLE '{key}'")

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self.list())

    def __len__(self):
        """The len() command"""
        self.cursor.execute(
            "SELECT COUNT(*)"
            "  FROM sqlite_master"
            " WHERE type = 'table'"
        )
        return self.cursor.fetchone()[0]

    def __repr__(self):
        """The string representation of this object"""
        raise NotImplementedError()

    def __str__(self):
        """The pretty string representation of this object"""
        raise NotImplementedError()

    def __contains__(self, item):
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.

        tmp_item = item.strip('"')
        self.cursor.execute(
            "SELECT COUNT(*)"
            "  FROM sqlite_master"
            f" WHERE type = 'table' and name = '{tmp_item}'"
        )
        return self.cursor.fetchone()[0] == 1

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    def list(self):
        """Return a list of all the tables in the system."""
        result = []
        for row in self.db.execute(
            "SELECT name"
            "  FROM sqlite_master"
            " WHERE type = 'table'"
        ):
            result.append(row['name'])
        return result

    @property
    def filename(self):
        """The name of the file (or URI) for the database."""
        return self._filename

    @filename.setter
    def filename(self, value):
        if value != self._filename:
            if self._db is not None:
                self.cursor.close()
                self._db.commit()
                self._db.close()
                self._db = None
                self._cursor = None
            self._filename = value
            if self._filename is not None:
                self._db = sqlite3.connect(self._filename)
                self._db.row_factory = sqlite3.Row
                self._db.execute('PRAGMA foreign_keys = ON')
                self._cursor = self._db.cursor()
                self._initialize()

    @property
    def db(self):
        return self._db

    @property
    def cursor(self):
        return self._cursor

    @property
    def n_atoms(self):
        """The number of atoms in the system."""
        return self.atoms.n_atoms

    @property
    def atoms(self):
        """The atoms, which are held as a dictionary of arrays"""
        return self['atoms']

    @property
    def n_bonds(self):
        """The number of bonds in the system."""
        return self.bonds.n_bonds

    @property
    def bonds(self):
        """The bonds, which are held as a dictionary of arrays"""
        return self['bonds']

    @property
    def cell(self):
        """The periodic cell."""
        return self['cell']

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        self.cursor.execute(
            "SELECT value FROM metadata WHERE item = 'periodicity'"
        )
        return int(self.cursor.fetchone()[0])

    @periodicity.setter
    def periodicity(self, value):
        if value < 0 or value > 3:
            raise ValueError('The periodicity must be between 0 and 3.')
        self.cursor.execute(
            "UPDATE metadata SET value = ? WHERE item = 'periodicity'",
            (str(value))
        )

    @property
    def version(self):
        """The version of the system, incrementing from 0"""
        self.cursor.execute(
            "SELECT value FROM metadata WHERE item = 'version'"
        )
        return int(self.cursor.fetchone()[0])

    @version.setter
    def version(self, value):
        self.cursor.execute(
            "UPDATE metadata SET value = ? WHERE item = 'version'",
            (str(value),)
        )

    @property
    def coordinate_type(self):
        """The coordinates of the system, 'fractional' or 'Cartesian'"""
        self.cursor.execute(
            "SELECT value FROM metadata WHERE item = 'coordinate type'"
        )
        return self.cursor.fetchone()[0]

    @coordinate_type.setter
    def coordinate_type(self, value):
        if value.lower()[0] == 'f':
            self.cursor.execute(
                "UPDATE metadata SET text = 'fractional'"
                " WHERE item = 'coordinate type'"
            )
        else:
            self.cursor.execute(
                "UPDATE metadata SET text = 'Cartesian'"
                " WHERE item = 'coordinate type'"
            )

    def create_table(self, name, cls=_Table, other=None):
        """Create a new table with the given name.

        Parameters
        ----------
        name : str
            The name of the new table.

        cls : _Table subclass
            The class of the new table, defaults to _Table

        Returns
        -------
        table : class Table
            The new table
        """
        if name in self:
            raise KeyError(f"'{name}' already exists in the system.")

        return cls(self, name, other)

    def _log_changes(self, previous, backup_file):
        """Track changes to the system"""
        changed = False

        # Attach the previous database in order to do comparisons.
        self.cursor.execute(f"ATTACH DATABASE '{backup_file}' AS previous")

        # Tables
        self.cursor.execute(
            "SELECT name"
            "  FROM sqlite_master"
            " WHERE type = 'table'"
        )
        tables = [x[0] for x in self.cursor.fetchall()]

        self.cursor.execute(
            "SELECT name"
            "  FROM previous.sqlite_master"
            " WHERE type = 'table'"
        )
        previous_tables = [x[0] for x in self.cursor.fetchall()]

        for table in tables:
            if table not in previous_tables:
                changed = True
                print(f'Table {table} added to the system.')
        for table in previous_tables:
            if table not in tables:
                changed = True
                print(f'Table {table} deleted from the system.')

        # Now check each table
        for table in tables:
            attributes = self.attributes(table)
            columns = set(attributes)
            if table in previous_tables:
                previous_attributes = self.attributes(f'previous.{table}')
                previous_columns = set(previous_attributes)
                added = columns - previous_columns
                if len(added) > 0:
                    changed = True
                    print(
                        f'The following columns were added to table {table}:'
                    )
                    print('\t' + '\n\t'.join(added))
                removed = previous_columns - columns
                if len(removed) > 0:
                    changed = True
                    print(
                        'The following columns were removed from table '
                        f'{table}:'
                    )
                    print('\t' + '\n\t'.join(removed))
                in_common = previous_columns & columns
                if len(in_common) > 0:
                    column_def = 'rowid, ' + ', '.join(in_common)

                    # Diff the tables....
                    print(f"Diff'ing table = {table}")

                    self.cursor.execute(
                        f'SELECT COUNT(*) FROM previous.{table}'
                    )
                    n_rows = self.cursor.fetchone()[0]
                    print(f'   There were {n_rows} rows.')

                    self.cursor.execute(f'SELECT COUNT(*) FROM {table}')
                    n_rows = self.cursor.fetchone()[0]
                    print(f'   Now there are {n_rows} rows.')

                    self.cursor.execute(
                        f'SELECT {column_def} FROM'
                        '       ('
                        f'           SELECT {column_def} FROM previous.{table}'
                        '           EXCEPT'
                        f'           SELECT {column_def} FROM {table}'
                        '       )'
                        ' UNION ALL '
                        f'SELECT {column_def} FROM'
                        '       ('
                        f'           SELECT {column_def} FROM {table}'
                        '           EXCEPT'
                        f'           SELECT {column_def} FROM previous.{table}'
                        '       )'
                        ' ORDER BY rowid'
                    )
                    lines = self.cursor.fetchall()
                    print(f'Found {len(lines)} lines of differences.')
                    if len(lines) > 0:
                        changed = True
                    last = None
                    if table == 'metadata':
                        for line in lines:
                            if last is None:
                                last = line
                            elif line['rowid'] == last['rowid']:
                                for k1, v1, v2 in zip(last.keys(), last, line):
                                    if v1 != v2:
                                        print(
                                            f"{line['item']} changed from "
                                            f"{v1} to {v2}"
                                        )
                            else:
                                last = line
                    else:
                        for line in lines:
                            if last is None:
                                last = line
                            elif line['rowid'] == last['rowid']:
                                for k1, v1, v2 in zip(last.keys(), last, line):
                                    if v1 != v2:
                                        print(
                                            f"{line['rowid']}, {k1} changed "
                                            f"from {v1} to {v2}"
                                        )
                            else:
                                last = line
        # Detach the previous database
        self.cursor.execute("DETACH DATABASE previous")

        if changed:
            self.version += 1
            print('The system was changed')
        else:
            print('The system was not changed')

        return changed

    def _initialize(self):
        """Initialize the SQLite database."""
        if 'elements' not in self:
            self._initialize_elements()

        if 'metadata' not in self:
            self._initialize_metadata()

    def _initialize_metadata(self):
        """Set up the table of metadata."""
        table = self['metadata']
        table.add_attribute('item', coltype='str')
        table.add_attribute('value', coltype='str')

        table.append(item='version', value='0')
        table.append(item='periodicity', value='0')
        table.append(item='coordinate type', value='Cartesian')

        self.db.commit()

    def _initialize_elements(self):
        """Set up the table of elements."""
        table = self['elements']
        table.add_attribute('atno', coltype='int', pk=True)
        table.add_attribute('symbol', coltype='str', index='unique')
        table.add_attribute('mass', coltype='float')

        for symbol, data in seamm_util.element_data.items():
            table.append(
                symbol=symbol,
                atno=data['atomic number'],
                mass=data['atomic weight']
            )
        self.db.commit()

    def attributes(self, tablename: str) -> Dict[str, Any]:
        """The attributes -- columns -- of a given table.

        Parameters
        ----------
        tablename : str
            The name of the table, optionally including the schema followed by
            a dot.

        Returns
        -------
        attributes : Dict[str, Any]
            A dictionary of dictionaries for the attributes and their
            descriptors
        """
        if '.' in tablename:
            schema, tablename = tablename.split('.')
            sql = f"PRAGMA {schema}.table_info('{tablename}')"
        else:
            sql = f"PRAGMA table_info('{tablename}')"

        result = {}
        for line in self.db.execute(sql):
            result[line['name']] = {
                'type': line['type'],
                'notnull': bool(line['notnull']),
                'default': line['dflt_value'],
                'primary key': bool(line['pk'])
            }
        return result


if __name__ == '__main__':  # pragma: no cover
    system = System()
    with system as s:
        s.periodicity = 3
        # s.coordinate_type = 'f'

    print(f'metadata? {"metadata" in system}')
    print(f'  table1? {"table1" in system}')

    metadata = system['metadata']
    # pprint.pprint(metadata.attributes)

    table1 = system.create_table('table1')
    table1.add_attribute('atno', coltype='int', default=-1)
    print('tables: ' + ', '.join(iter(system)))

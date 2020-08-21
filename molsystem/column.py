# -*- coding: utf-8 -*-

import collections.abc
import itertools
import logging
# import pandas
from typing import Any

# from molsystem.system import System
# from molsystem.table import _Table

# tpSystem = TypeVar("tpSystem", System, None)
# tpTable = TypeVar("tpTable", _Table, None)

logger = logging.getLogger(__name__)


class _Column(collections.abc.MutableSequence):
    """A list-like object for holding a column data

    This is a wrapper around a single colmun in a SQL table. Setting a value
    updates the SQL database appropriately.
    """

    def __init__(self, table, column: str) -> None:
        self._table = table
        self._column = column
        self._data = None
        self._rowids = None
        self._initialize()

    def __getitem__(self, index) -> Any:
        """Allow [] to access the data!"""
        print(index)
        print(type(index))
        return self._data[index]

    def __setitem__(self, index, value) -> None:
        """Allow x[index] access to the data"""
        length = self._table.length_of_values(value)
        db = self._table.db
        table = self._table.table
        if isinstance(index, int):
            if length == 0:
                # Scalar
                db.execute(
                    f"UPDATE {table} SET {self.column} = ? WHERE rowid = ?",
                    (value, self._rowids[index])
                )
            elif length == 1:
                # One-element list
                db.execute(
                    f"UPDATE {table} SET {self.column} = ? WHERE rowid = ?",
                    (value[0], self._rowids[index])
                )
            else:
                raise IndexError('Only 1 value required to update 1 item')
        elif isinstance(index, slice):
            rowids = []
            for rowid in itertools.islice(
                self._rowids, index.start, index.stop, index.step
            ):
                rowids.append(rowid)
            parameters = []
            if length == 0:
                for rowid in rowids:
                    parameters.append((value, rowid))
            elif length == 1:
                for rowid in rowids:
                    parameters.append((value[0], rowid))
            elif length == len(rowids):
                for rowid, val in zip(rowids, value):
                    parameters.append((val, rowid))
            else:
                raise IndexError(
                    f'The number of values ({length}) must be 1 or the size '
                    f'of the slice ({len(rowids)}).'
                )
            db.executemany(
                f"UPDATE {table} SET {self.column} = ? WHERE rowid = ?",
                parameters
            )

    def __delitem__(self, index, value) -> None:
        """Do NOT allow deletion!"""
        raise RuntimeError('Items may not be removed from a Column')

    def __len__(self) -> int:
        """The len() command"""
        return len(self._data)

    def insert(self, index, value) -> None:
        """Do NOT allow insertion!"""
        raise RuntimeError('Items may not be inserted into a Column')

    def __eq__(self, other):
        """If this is equal to another list."""
        print("in column __eq__")
        print(self._data)
        print(other)
        return self._data == other

    @property
    def column(self):
        """A nicely quoted version of the column name"""
        return '"' + self._column + '"'

    def _initialize(self):
        """Get the column data from the underlying table and cache it."""
        self._data = []
        self._rowids = []
        print(f'SELECT rowid, {self.column} FROM {self._table.table}')
        for row in self._table.db.execute(
            f'SELECT rowid, {self.column} FROM {self._table.table}'
        ):
            print(row)
            for i in row:
                print(i)
            self._rowids.append(row[0])
            self._data.append(row[1])
        print(self._rowids)
        print(self._data)


if __name__ == '__main__':  # pragma: no cover
    import os.path
    import tempfile

    from molsystem import System

    x = [1.0, 2.0, 3.0]
    y = [4.0, 5.0, 6.0]
    z = [7.0, 8.0, 9.0]
    atno = [8, 1, 1]

    with tempfile.TemporaryDirectory() as tmpdirname:
        print(tmpdirname)
        filepath = os.path.join(tmpdirname, 'system1.db')
        system1 = System(filename=filepath)

        table1 = system1.create_table('table1')
        table1.add_attribute('atno', coltype='int')
        for column in ('x', 'y', 'z'):
            table1.add_attribute(column, coltype='float')
        with table1 as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno)

        atno = table1['atno']
        print(atno[1])
        print(atno[0:2])

# -*- coding: utf-8 -*-

import collections.abc
import itertools
import logging

from molsystem.frozencolumn import _FrozenColumn as FrozenColumn

logger = logging.getLogger(__name__)


class _Column(FrozenColumn, collections.abc.MutableSequence):
    """A list-like object for holding a column data

    This is a wrapper around a single column in a SQL table. Setting a value
    updates the SQL database appropriately.
    """

    def __init__(self, table, column: str, sql=None, where=""):
        self._rowids = None
        super().__init__(table, column, sql, where)

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
                    (value, self._rowids[index]),
                )
                self._data[index] = value
            elif length == 1:
                # One-element list
                db.execute(
                    f"UPDATE {table} SET {self.column} = ? WHERE rowid = ?",
                    (value[0], self._rowids[index]),
                )
                self._data[index] = value[0]
            else:
                raise IndexError("Only 1 value required to update 1 item")
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
                self._data[index] = [value] * len(rowids)
            elif length == 1:
                for rowid in rowids:
                    parameters.append((value[0], rowid))
                self._data[index] = [value[0]] * len(rowids)
            elif length == len(rowids):
                for rowid, val in zip(rowids, value):
                    parameters.append((val, rowid))
                self._data[index] = value
            else:
                raise IndexError(
                    f"The number of values ({length}) must be 1 or the size "
                    f"of the slice ({len(rowids)})."
                )
            db.executemany(
                f"UPDATE {table} SET {self.column} = ? WHERE rowid = ?", parameters
            )
        db.commit()

    def __delitem__(self, index, value) -> None:
        """Do NOT allow deletion!"""
        raise RuntimeError("Items may not be removed from a Column")

    def __repr__(self) -> str:
        """The string representation of this object"""
        return repr([(r, x) for r, x in zip(self._rowids, self._data)])

    def insert(self, index, value) -> None:
        """Do NOT allow insertion!"""
        raise RuntimeError("Items may not be inserted into a Column")

    def _initialize(self):
        """Get the column data from the underlying table and cache it."""
        self._data = []
        self._rowids = []
        if self._sql is None:
            sql = (
                f"SELECT rowid, {self.column}"
                f"  FROM {self._table.table} {self._where}"
            )
        else:
            sql = self._sql

        for row in self._table.db.execute(sql):
            self._rowids.append(row[0])
            self._data.append(row[1])

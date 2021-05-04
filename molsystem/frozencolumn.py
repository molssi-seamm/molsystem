# -*- coding: utf-8 -*-

import collections.abc
import logging

logger = logging.getLogger(__name__)


class _FrozenColumn(collections.abc.Sequence):
    """A list-like object for holding a column data

    This is a wrapper around a single column in a SQL table.
    """

    def __init__(self, table, column: str, sql=None, where=""):
        self._table = table
        self._column = column
        self._sql = sql
        self._where = where
        self._data = None

        self._initialize()

    def __getitem__(self, index):
        """Allow [] to access the data!"""
        return self._data[index]

    def __len__(self) -> int:
        """The len() command"""
        return len(self._data)

    def __repr__(self) -> str:
        """The string representation of this object"""
        return repr(self._data)

    def __str__(self) -> str:
        """The pretty string representation of this object"""
        return str(self._data)

    def __eq__(self, other):
        """If this is equal to another list."""
        return self._data == other

    @property
    def column(self):
        """A nicely quoted version of the column name"""
        return '"' + self._column + '"'

    def _initialize(self):
        """Get the column data from the underlying table and cache it."""
        self._data = []
        if self._sql is None:
            sql = f"SELECT {self.column} FROM {self._table.table} {self._where}"
        else:
            sql = sql._sql

        for row in self._table.db.execute(sql):
            self._data.append(row[0])

    def equal(self, other, tol=1.0e-6):
        """Check if we are equal to another iterable to within a tolerance.

        Parameters
        ----------
        other : iterable
            The other object to check against
        tol : float = 1.0e-06
            The tolerance for comparing floating point numbers.

        Returns
        -------
        equals : bool
            Boolean indicating whether the two are equal.
        """
        if len(self._data) != len(other):
            return False

        if len(self._data) == 0:
            return True

        if isinstance(self._data[0], float):
            for v1, v2 in zip(self._data, other):
                if abs(v1 - v2) > tol:
                    return False
            return True
        else:
            return self._data == other

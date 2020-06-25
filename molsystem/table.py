# -*- coding: utf-8 -*-

import collections.abc
import copy
import logging
from typing import Any, Dict, List

import numpy as np
import pandas as pd
logger = logging.getLogger(__name__)


class Table(collections.abc.MutableMapping):
    """A dictionary-like object for holding tabular data

    For efficiency the data is stored as arrays using Pandas. In addtion, there
    may be empty ("free") space at the end of the table to make adding small
    numbers of rows repeatedly reasonable efficient.

    Metadata is used and stored to define the attributes ("columns").

    Attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Rows can be added ('append') or removed ('delete').
    """

    allocate_min = 128
    """The minimum number of rows to allocate at a time."""

    allocate_factor = 0.5
    """Allocate the current number of rows times this factor -- but at least
    'allocate_min' rows -- when space is needed."""

    def __init__(self, other: 'Table' = None) -> None:

        # this is the only thing not copied
        self._checkpoints = []
        self._word = 'rows'

        # copy constructor
        if other and isinstance(other, Table):
            self._public = copy.deepcopy(other._public)
            self._private = copy.deepcopy(other._private)
        else:
            self._public = pd.DataFrame(
                index=pd
                .RangeIndex(start=0, stop=Table.allocate_min, name='uid')
            )
            self._private = {
                'version': 0,
                'n_rows': 0,
                'attributes': {},
            }

    def __enter__(self) -> Any:
        tmp = copy.deepcopy(self)
        self._checkpoints.append(tmp)
        return tmp

    def __exit__(self, etype, value, traceback) -> None:
        if etype is None:
            # No exception occurred, so replace ourselves with the tmp copy
            tmp = self._checkpoints.pop()
            self._public, tmp._public = tmp._public, self._public
            self._private, tmp._private = tmp._private, self._private

            # and log the changes
            self._log_changes(tmp)

    def __getitem__(self, key) -> Any:
        """Allow [] access to the dictionary!"""
        return self._public[key][0:self.n_rows]

    def __setitem__(self, key, value) -> None:
        """Allow x[key] access to the data"""
        self._public[key][0:self.n_rows] = value

    def __delitem__(self, key) -> None:
        """Allow deletion of keys"""
        self._public.drop(key, axis=1, inplace=True)

    def __iter__(self) -> Any:
        """Allow iteration over the object"""
        return iter(self._public)

    def __len__(self) -> int:
        """The len() command"""
        return len(self._public)

    def __repr__(self) -> str:
        """The string representation of this object"""
        return repr(self._public[0:self.n_rows])

    def __str__(self) -> str:
        """The pretty string representation of this object"""
        return self._public[0:self.n_rows].to_string()

    def __contains__(self, item) -> bool:
        """Return a boolean indicating if a key exists."""
        if item in self._public:
            return True
        return False

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        return self._public.equals(other._public)

    @property
    def version(self) -> int:
        """The version of the table, which increments monotonically as changes
        are made."""
        return self._private['version']

    @property
    def n_rows(self) -> int:
        """The number of rows in the table."""
        return self._private['n_rows']

    @property
    def attributes(self) -> Dict[str, Any]:
        """The definitions of the attributes."""
        return self._private['attributes']

    @property
    def free(self) -> int:
        """How many rows are pre-allocated but free.

        When setting this, space will be allocated if needed so that this
        number of rows can be added without further allocation.
        """

        return len(self) - self.n_rows

    @free.setter
    def free(self, n: int) -> None:
        n_add = n - self.free
        if n_add > 0:
            self.allocate(n_add)

    @property
    def index(self) -> pd.RangeIndex:
        """The index of the underlying pandas table."""
        return self._public.index

    def _log_changes(self, previous: pd.DataFrame) -> bool:
        """Track changes to the table"""
        changed = False
        self._private['version'] += 1
        for key in self:
            if key not in previous:
                print("attribute '{}' added".format(key))
                changed = True
            elif not np.array_equal(self[key], previous[key]):
                print("'{}' changed.".format(key))
                changed = True
        for key in previous:
            if key not in self:
                print("attribute '{}' deleted".format(key))
                changed = True

        if not changed:
            self._private['version'] -= 1
            print('The table was not changed')

        return changed

    def add_attribute(
        self,
        name: str,
        coltype: Any = None,
        default: Any = None,
        values: List = None
    ) -> None:
        """Adds a new attribute.

        If the attribute has been defined previously, the column type and
        default will automatically be those given in the definition. If they
        are provided as arguments, they must be the same as in the definition.

        If the default value is None, you must always provide values wherever
        needed, for example when adding a row.

        Args:
            name: the name of the attribute.
            coltype: the type of the attribute (column). Must be an acceptable
                numpy type.
            default: the default value for the attribute if no value is given.
            values: either a single value or a list of values length 'nrows' of
                values to fill the column.

        Returns:
            None
        """

        if name in self:
            raise RuntimeError(
                "{} attribute '{}' is already defined!".format(
                    self.__class__.__name__, name
                )
            )

        if name not in self.attributes:
            # Not yet defined
            attr = self.define_attribute(name, coltype, default)
        else:
            # If we have the attributes, check that any given are reasonable.
            attr = self.attributes[name]
            if coltype is not None:
                if attr['type'] != coltype:
                    raise ValueError(
                        "Column type should be '{}', not '{}'".format(
                            attr['type'], coltype
                        )
                    )
            if default is not None:
                if default != attr['default']:
                    raise ValueError(
                        "Default should be '{}', not '{}'".format(
                            attr['default'], default
                        )
                    )

        n = len(self)
        # see if the values are given
        if values is not None:
            try:
                length = len(values)
            except TypeError:
                length = 1

            if length != 1 and length != self.n_rows:
                raise IndexError(
                    "The number of values given, "
                    "{}, must be either 1, or the number of {}: {}".format(
                        length, self._word, self.n_rows
                    )
                )

        column = np.full(
            shape=(n,), fill_value=attr['default'], dtype=attr['type']
        )

        if values is not None:
            column[0:self.n_rows] = values

        added = pd.DataFrame(data=column, columns=[name])
        self._public = self._public.join(added)

    def define_attribute(
        self,
        name: str,
        coltype: Any = None,
        default: Any = None
    ) -> Dict[str, Any]:
        """Defines a new attribute

        If the default value is None, you must always provide values wherever
        needed, for example when adding a row.

        Args:
            name: the name of the attribute.
            coltype: the type of the attribute (column). Must be an acceptable
                numpy type.
            default: the default value for the attribute if no value is given.

        Returns:
            None
        """

        if name in self.attributes:
            raise RuntimeError(
                "{} attribute '{}' is already defined!".format(
                    self.__class__.__name__, name
                )
            )

        # Setup the attributes
        tmp = self.attributes[name] = {}
        tmp['type'] = coltype
        tmp['default'] = default

        return tmp

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more rows

        The keywords are the names of attributes and the value to use.
        The default value for any attributes not given is used unless it is
        'None' in which case an error is thrown. It is an error if there is not
        an exisiting attribute corresponding to any given as arguments.

        Args:
            kwargs: any number <attribute name> = <value> keyword arguments
                giving existing attributes and values.

        Returns:
            None
        """

        # Check keys and lengths of added rows
        n_rows = None

        for key, value in kwargs.items():
            if key not in self:
                raise KeyError(
                    '"{key}" is not an attribute of the {cls}!'.format(
                        key=key, cls=self.__class__.__name__.lower()
                    )
                )
            try:
                length = len(value)
            except TypeError:
                length = 1

            if n_rows is None:
                n_rows = length

            if length != 1 and length != n_rows:
                if n_rows == 1:
                    n_rows = length
                else:
                    raise IndexError(
                        'key "{}" has the wrong number of values, '
                        .format(key) +
                        '{}. Should be 1 or the number of {} ({}).'
                        .format(length, self._word, n_rows)
                    )

        # Check that any missing attributes have defaults
        for key in self:
            if key not in kwargs:
                if key not in self.attributes:
                    # Don't think this can happen ... but test anyway
                    raise KeyError(
                        "Cannot handle attribute '{}'.".format(key)
                    )  # pragma: no cover
                if 'default' not in self.attributes[key] or \
                   self.attributes[key]['default'] is None:
                    raise KeyError(
                        "There is no default for attribute "
                        "'{}'. You must supply a value".format(key)
                    )

        # All okay, so proceed.
        to_add = {}
        for key, value in kwargs.items():
            to_add[key] = np.full(
                shape=(n_rows,),
                fill_value=self.attributes[key]['default'],
                dtype=self.attributes[key]['type']
            )
            to_add[key] = value
        for key in self:
            if key not in kwargs:
                to_add[key] = np.full(
                    shape=(n_rows,),
                    fill_value=self.attributes[key]['default'],
                    dtype=self.attributes[key]['type']
                )

        df_tmp = pd.DataFrame(
            to_add,
            index=pd.RangeIndex(
                start=self.n_rows, stop=self.n_rows + n_rows, name='uid'
            )
        )

        # check that there is space
        if n_rows > self.free:
            delta = (
                int(Table.allocate_factor * (self.n_rows + n_rows)) - self.free
            )
            if delta < Table.allocate_min:
                delta = Table.allocate_min
            self.allocate(delta)

        self._public.update(df_tmp)
        self._private['n_rows'] += n_rows

    def allocate(self, n: int) -> None:
        """Pre-allocates space for <n> more rows.

        Ensures that there is space for <n> more rows. Note that if there is
        already sufficient space, nothing is done. This method will not reduce
        the available space. See 'trim' for reducing space.

        Args:
            n: the number of rows to add, >= 1

        Returns:
            None
        """
        logger.debug('Allocating space for {} more {}'.format(n, self._word))

        to_add = {}
        for key in self:
            default = self.attributes[key]['default']
            if default == 0:
                to_add[key] = np.zeros(
                    shape=(n,), dtype=self.attributes[key]['type']
                )
            else:
                to_add[key] = np.full(
                    shape=(n,),
                    fill_value=default,
                    dtype=self.attributes[key]['type']
                )

        start = self.index.max() + 1
        df_tmp = pd.DataFrame(
            to_add,
            index=pd.RangeIndex(start=start, stop=start + n, name='uid')
        )
        self._public = self._public.append(df_tmp)

    def trim(self) -> None:
        """Removes any extra space at the end of the table.

        When rows are added, space is allocated for a reasonable number of rows
        for performance reasons. The method removes any unused space, trimming
        the table to exactly the needed size.

        Args:
            None

        Returns:
            None
        """
        if self.free > 0:
            self._public.drop(self.index.to_list()[self.n_rows:], inplace=True)


if __name__ == '__main__':  # pragma: no cover
    import numpy.random as nprand
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False) -> None:
        table = Table()
        for column in ('a', 'b', 'c'):
            table.add_attribute(column, coltype=np.float64, default=np.nan)

        a = nprand.uniform(low=0, high=100, size=nper)
        b = nprand.uniform(low=0, high=100, size=nper)
        c = nprand.uniform(low=0, high=100, size=nper)

        with table as tmp:
            if preallocate:
                tmp.allocate(nper * nrepeat)
            for i in range(0, nrepeat):
                tmp.append(a=a, b=b, c=c)
                a += 10.0
                b += 5.0
                c -= 3.0

    nrepeat = 1000
    nper = 100
    t = timeit.timeit(
        "run(nper={}, nrepeat={})".format(nper, nrepeat),
        setup="from __main__ import run",
        timer=time.time,
        number=1
    )

    print("Creating {} rows took {:.3f} s".format(nper * nrepeat, t))

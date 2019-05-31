# -*- coding: utf-8 -*-

import collections.abc
import copy
import logging
import numpy as np
import pandas as pd
import pprint

"""A dictionary-like object for holding a system
"""

logger = logging.getLogger(__name__)


class System(collections.abc.MutableMapping):
    def __init__(self, **kwargs):
        self._checkpoints = []
        self._version = 0

        self._data = dict()

        self['atoms'] = pd.DataFrame(
            {
                'x': np.zeros((0,)),
                'y': np.zeros((0,)),
                'z': np.zeros((0,)),
                'atno': np.zeros((0,), dtype=np.int8)
            }
        )
        self['periodicity'] = 0
        self['coordinates'] = 'Cartesian'

        self._data.update(**kwargs)

    def __enter__(self):
        tmp = copy.deepcopy(self)
        self._checkpoints.append(tmp)
        return tmp

    def __exit__(self, etype, value, traceback):
        if etype is None:
            # No exception occurred, so log changes
            tmp = self._checkpoints.pop()
            self._log_changes(tmp)

    def __getitem__(self, key):
        """Allow [] access to the dictionary!"""
        return self._data[key]

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        self._data[key] = value

    def __delitem__(self, key):
        """Allow deletion of keys"""
        del self._data[key]

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self._data)

    def __len__(self):
        """The len() command"""
        return len(self._data)

    def __repr__(self):
        """The string representation of this object"""
        return repr(self._data)

    def __str__(self):
        """The pretty string representation of this object"""
        return pprint.pformat(self._data)

    def __contains__(self, item):
        """Return a boolean indicating if a key exists."""
        if item in self._data:
            return True
        return False

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        return self._data == other._data

    def copy(self):
        """Return a shallow copy of the dictionary"""
        return self._data.copy()

    @property
    def version(self):
        """The version of the system, which increments monotonically"""
        return self._version

    @property
    def atoms(self):
        """The atoms, which are held as a dictionary of arrays"""
        return self._data['atoms']

    @atoms.setter
    def atoms(self, atoms):
        self._data['atoms'] = atoms

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        return self._data['periodicity']

    @periodicity.setter
    def periodicity(self, value):
        try:
            if value < 0 or value > 3:
                raise ValueError(
                    "The periodicity must be 0, 1, 2 or 3, not '{}'"
                    .format(value)
                )
        except TypeError:
            raise ValueError(
                "The periodicity must be an integer 0, 1, 2 or 3, not '{}'"
                .format(value)
            )
        self._data['periodicity'] = value

    @property
    def coordinates(self):
        """The coordinates of the system, 'fractional' or 'Cartesian'"""
        return self._data['coordinates']

    @coordinates.setter
    def coordinates(self, value):
        try:
            if 'fractional'.startswith(value.lower()):
                self._data['coordinates'] = 'fractional'
            elif 'cartesian'.startswith(value.lower()):
                self._data['coordinates'] = 'Cartesian'
            else:
                raise ValueError(
                    "The coordinates must be 'Cartesian' or 'fractional', "
                    "not '{}'".format(value)
                )
        except TypeError:
            raise ValueError(
                "The coordinates must be 'Cartesian' or 'fractional', "
                "not '{}'".format(value)
            )

    def _log_changes(self, other):
        """Track changes to the system"""
        changed = False
        self._version += 1
        for key in self:
            if key == 'atoms':
                pass
            elif self[key] != other[key]:
                print("'{}' changed.".format(key))
                changed = True
        if not changed:
            self._version -= 1
            print('The system was not changed')


if __name__ == '__main__':
    system = System()
    with system as s:
        s.periodicity = 3
        s.coordinates = 'f'

# -*- coding: utf-8 -*-

"""A dictionary-like object for holding atoms

For efficiency the atom-data is stored as arrays -- numpy if possible, python
lists otherwise, along with metadata to define the attributes ("columns").

In some ways this a bit like pandas; however, we also need more control than a
simple pandas dataframe provides....
"""

import copy
import logging
from typing import Any, Dict, TypeVar

import numpy as np

import molsystem

System_tp = TypeVar("System_tp", "System", "Atoms", None)
Atoms_tp = TypeVar("Atoms_tp", "Atoms", None)

logger = logging.getLogger(__name__)

attributes = {
    'x': {
        'type': np.float64,
        'default': 0.0,
    },
    'y': {
        'type': np.float64,
        'default': 0.0,
    },
    'z': {
        'type': np.float64,
        'default': 0.0,
    },
    'atno': {
        'type': np.int8,
        'default': -1,
    },
    'name': {
        'type': np.object,
        'default': '',
    },
}


class Atoms(molsystem.Table):
    """The Atoms class holds arrays of attributes describing atoms

    Two attributes are required: the coordinates ('xyz') and atomic numbers
    ('atno'). Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Atoms can be added ('append') or removed ('delete').
    """

    def __init__(self, other: System_tp = None) -> None:

        super().__init__()

        self._word = 'atoms'

        # copy constructor
        if isinstance(other, Atoms):
            self._coordinate_type = other._coordinate_type
            self._public = copy.deepcopy(other._public)
            self._private = copy.deepcopy(other._private)
        else:
            self._coordinate_type = 'Cartesian'
            self._private['system'] = other
            self._private['attributes'] = attributes
            for key in ('x', 'y', 'z', 'atno'):
                self.add_attribute(key)

    def __exit__(self, etype, value, traceback) -> None:
        if etype is None:
            # No exception occurred, so replace ourselves with the tmp copy
            tmp = self._checkpoints.pop()
            self._coordinate_type, tmp._coordinate_type = (
                tmp._coordinate_type, self._coordinate_type
            )
            self._public, tmp._public = tmp._public, self._public
            self._private, tmp._private = tmp._private, self._private

            # and log the changes
            self._log_changes(tmp)

    @property
    def n_atoms(self) -> int:
        """The number of atoms"""
        return self._private['n_rows']

    @property
    def coordinate_type(self):
        """The type of coordinates: 'fractional' or 'Cartesian'"""
        return self._coordinate_type

    @coordinate_type.setter
    def coordinate_type(self, value):
        try:
            if 'fractional'.startswith(value.lower()):
                self._coordinate_type = 'fractional'
            elif 'cartesian'.startswith(value.lower()):
                self._coordinate_type = 'Cartesian'
            else:
                raise ValueError(
                    "The coordinate_type must be 'Cartesian' or 'fractional', "
                    "not '{}'".format(value)
                )
        except Exception:
            raise ValueError(
                "The coordinate_type must be 'Cartesian' or 'fractional', "
                "not '{}'".format(value)
            )

    def _log_changes(self, previous: Atoms_tp) -> bool:
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
        if self._coordinate_type != previous._coordinate_type:
            changed = True
            print('type of coordinates changed')

        if not changed:
            self._private['version'] -= 1
            print('The Atoms object not changed')

        return changed

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more atoms

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Check keys and lengths of added atoms
        if 'x' not in kwargs or 'y' not in kwargs or 'z' not in kwargs:
            raise KeyError("The coordinates are required!")

        super().append(**kwargs)


if __name__ == '__main__':  # pragma: no cover
    import numpy.random as nprand
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False) -> None:
        system = None
        atoms = Atoms(system)
        x = nprand.uniform(low=0, high=100, size=nper)
        y = nprand.uniform(low=0, high=100, size=nper)
        z = nprand.uniform(low=0, high=100, size=nper)
        atno = np.array(nper * [6])
        with atoms as tmp:
            if preallocate:
                tmp.allocate(nper * nrepeat)
            for i in range(0, nrepeat):
                tmp.append(x=x, y=y, z=z, atno=atno)
                x += 10.0
                y += 5.0
                z -= 3.0

    nrepeat = 1000
    nper = 100
    t = timeit.timeit(
        "run(nper={}, nrepeat={})".format(nper, nrepeat),
        setup="from __main__ import run",
        timer=time.time,
        number=1
    )

    print("Creating {} atoms took {:.3f} s".format(nper * nrepeat, t))

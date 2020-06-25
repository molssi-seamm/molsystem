# -*- coding: utf-8 -*-

import collections.abc
import copy
import logging
import pprint
from typing import Any

import molsystem
"""A dictionary-like object for holding a system
"""

logger = logging.getLogger(__name__)


class System(collections.abc.MutableMapping):

    def __init__(self, other=None, **kwargs):
        # this is the only thing not copied
        self._checkpoints = []

        # copy constructor
        if other and isinstance(other, System):
            self._public = copy.deepcopy(other._public)
            self._private = copy.deepcopy(other._private)
        else:
            self._public = {}
            self._public['atoms'] = molsystem.Atoms()
            self._public['bonds'] = molsystem.Bonds()
            self._public['cell'] = molsystem.Cell()

            self._private = {
                'version': 0,
            }
        self._public.update(**kwargs)

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

            # and delete the copy
            del tmp

    def __getitem__(self, key):
        """Allow [] access to the dictionary!"""
        return self._public[key]

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        self._public[key] = value

    def __delitem__(self, key):
        """Allow deletion of keys"""
        del self._public[key]

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self._public)

    def __len__(self):
        """The len() command"""
        return len(self._public)

    def __repr__(self):
        """The string representation of this object"""
        return repr(self._public)

    def __str__(self):
        """The pretty string representation of this object"""
        return pprint.pformat(self._public)

    def __contains__(self, item):
        """Return a boolean indicating if a key exists."""
        return item in self._public

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        return self._public == other._public

    @property
    def version(self):
        """The version of the system, which increments monotonically"""
        return self._private['version']

    @property
    def n_atoms(self):
        """The number of atoms in the system."""
        return self.atoms.n_atoms

    @property
    def atoms(self):
        """The atoms, which are held as a dictionary of arrays"""
        return self._public['atoms']

    @property
    def n_bonds(self):
        """The number of bonds in the system."""
        return self.bonds.n_bonds

    @property
    def bonds(self):
        """The bonds, which are held as a dictionary of arrays"""
        return self._public['bonds']

    @property
    def cell(self):
        """The periodic cell."""
        return self._public['cell']

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        return self['cell'].periodicity

    @periodicity.setter
    def periodicity(self, value):
        self['cell'].periodicity = value

    @property
    def coordinate_type(self):
        """The coordinates of the system, 'fractional' or 'Cartesian'"""
        return self.atoms.coordinate_type

    @coordinate_type.setter
    def coordinate_type(self, value):
        self.atoms.coordinate_type = value

    def _log_changes(self, previous):
        """Track changes to the system"""
        changed = False
        self._private['version'] += 1
        for key in self:
            if key in ('atoms', 'bonds', 'cell'):
                if self[key]._log_changes(previous[key]):
                    changed = True
            elif key not in previous:
                print(f"'{key}' added")
                changed = True
            elif self[key] != previous[key]:
                print(f"'{key}' changed.")
                changed = True
        for key in previous:
            if key not in self:
                print(f"'{key}' removed")
                changed = True

        if not changed:
            self._private['version'] -= 1
            print('The system was not changed')

        return changed


if __name__ == '__main__':  # pragma: no cover
    system = System()
    with system as s:
        s.periodicity = 3
        s.coordinate_type = 'f'

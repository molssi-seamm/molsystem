# -*- coding: utf-8 -*-

import collections.abc
import copy
import logging
from typing import Any, Iterable

logger = logging.getLogger(__name__)

labels = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')


class Cell(object):
    """The representation of the periodic cell
    """

    def __init__(self, other: 'Cell' = None) -> None:

        # this is the only thing not copied
        self._checkpoints = []

        # copy constructor
        if other and isinstance(other, Cell):
            self._public = copy.deepcopy(other._public)
            self._private = copy.deepcopy(other._private)
        else:
            self._public = {
                'periodicity': 0,
                'cell': [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
            }
            self._private = {
                'version': 0,
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

            # and delete the copy
            del tmp

    def __repr__(self) -> str:
        """The string representation of this object"""
        return repr(self._public)

    def __str__(self) -> str:
        """The pretty string representation of this object"""
        if self.periodicity == 0:
            return 'not periodic'
        elif self.periodicity == 3:
            return (
                f'3-D periodic, cell = {self.a:.3f} {self.b:.3f} {self.c:.3f} '
                f'{self.alpha:.1f} {self.beta:.1f} {self.gamma:.1f}'
            )

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        return self._public == other._public

    @property
    def version(self) -> int:
        """The version of the cell, which increments monotonically as changes
        are made."""
        return self._private['version']

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        return self._public['periodicity']

    @periodicity.setter
    def periodicity(self, value):
        if not isinstance(value, int):
            raise TypeError('The periodicity must be an integer.')

        if value < 0 or value > 3:
            raise ValueError(
                f"The periodicity must be 0, 1, 2 or 3, not '{value}'"
            )

        if value == 1 or value == 2:
            raise NotImplementedError(
                '1-D and 2-D periodicity not implemented yet.'
            )

        self._public['periodicity'] = value

    @property
    def parameters(self):
        """The list of cell parameters."""
        return self._public['cell']

    @parameters.setter
    def parameters(self, value: Iterable[float]):
        if not isinstance(value, collections.abc.Iterable):
            raise TypeError(
                'The cell parameters must be an iterable of floats, '
                'length 6.'
            )
        if len(value) != 6:
            raise TypeError(
                'The cell parameters must be an iterable of length 6.'
            )
        for val in value:
            if not isinstance(val, float) and not isinstance(val, int):
                raise TypeError('The cell parameters must be 6 floats.')

        self._public['cell'] = list(value)

    @property
    def a(self):
        """The length of the first lattice direction."""
        return self._public['cell'][0]

    @a.setter
    def a(self, value):
        self._public['cell'][0] = value

    @property
    def b(self):
        """The length of the second lattice direction."""
        return self._public['cell'][1]

    @b.setter
    def b(self, value):
        self._public['cell'][1] = value

    @property
    def c(self):
        """The length of the third lattice direction."""
        return self._public['cell'][2]

    @c.setter
    def c(self, value):
        self._public['cell'][2] = value

    @property
    def alpha(self):
        """The angle between b and c"""
        return self._public['cell'][3]

    @alpha.setter
    def alpha(self, value):
        self._public['cell'][3] = value

    @property
    def beta(self):
        """The angle between a and c"""
        return self._public['cell'][4]

    @beta.setter
    def beta(self, value):
        self._public['cell'][4] = value

    @property
    def gamma(self):
        """The angle between a and b"""
        return self._public['cell'][5]

    @gamma.setter
    def gamma(self, value):
        self._public['cell'][5] = value

    def _log_changes(self, previous: Any) -> bool:
        """Track changes to the table"""
        changed = False
        self._private['version'] += 1

        if self._public['periodicity'] != previous._public['periodicity']:
            changed = True
            print(f'periodicity changed to {self.periodicity}')
        for value, old, label in zip(
            self._public['cell'], previous._public['cell'], labels
        ):
            if value != old:
                print(f"cell parameter '{label}' changed")
                changed = True

        if not changed:
            self._private['version'] -= 1
            print('The cell was not changed')

        return changed

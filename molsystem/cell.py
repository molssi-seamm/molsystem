# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)


class Cell(object):
    """A class to handle cell parameters and their transformations."""

    def __init__(self, a, b, c, alpha, beta, gamma):
        self._parameters = [a, b, c, alpha, beta, gamma]

    def __getitem__(self, key):
        """Allow [] to access the data!"""
        return self._parameters[key]

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        self._parameters[key] = value

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self._parameters)

    def __len__(self) -> int:
        """The len() command"""
        return len(self._parameters)

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        return self.equal(other, tol=1.0e-12)

    def __repr__(self):
        """The representation of this object"""
        return repr(self._parameters)

    def __str__(self):
        """The pretty string representation of this object"""
        return str(self._parameters)

    @property
    def a(self):
        """The length of the first cell vector."""
        return self._parameters[0]

    @a.setter
    def a(self, value):
        self._parameters[0] = value

    @property
    def b(self):
        """The length of the second cell vector."""
        return self._parameters[1]

    @b.setter
    def b(self, value):
        self._parameters[1] = value

    @property
    def c(self):
        """The length of the third cell vector."""
        return self._parameters[2]

    @c.setter
    def c(self, value):
        self._parameters[2] = value

    @property
    def alpha(self):
        """The angle between b and c."""
        return self._parameters[3]

    @alpha.setter
    def alpha(self, value):
        self._parameters[3] = value

    @property
    def beta(self):
        """The angle between a and c."""
        return self._parameters[4]

    @beta.setter
    def beta(self, value):
        self._parameters[4] = value

    @property
    def gamma(self):
        """The angle between a and b."""
        return self._parameters[5]

    @gamma.setter
    def gamma(self, value):
        self._parameters[5] = value

    @property
    def parameters(self):
        """The cell parameters as a list."""
        return list(self._parameters)

    @parameters.setter
    def parameters(self, value):
        if len(value) != 6:
            raise ValueError('parameters must be of length 6')
        self._parameters = list(value)

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
        if len(other) != 6:
            return False

        for i, j in zip(self, other):
            if abs(i - j) > tol:
                return False

        return True

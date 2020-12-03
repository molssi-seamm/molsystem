# -*- coding: utf-8 -*-

import logging
import math

import numpy

logger = logging.getLogger(__name__)


def cos(value):
    return math.cos(math.radians(value))


def sin(value):
    return math.sin(math.radians(value))


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

    @property
    def volume(self):
        """The volume of the cell."""
        a, b, c, alpha, beta, gamma = self.parameters
        # Roundoff errors!
        value = cos(alpha) * cos(beta) * cos(gamma)
        if value < 0.0 and abs(value) < 1.0e-8:
            value = 0.0
        return (
            a * b * c * (1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2) +
            2 * math.sqrt(value)
        )

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

    def to_cartesians(self, uvw, as_array=False):
        """Convert fraction coordinates to Cartesians

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        uvw : [N][3*float] or ndarray
            The fractional coordinates.

        Returns
        -------
        xyz : [N][float*3] or ndarray
            The Cartesian coordinates.
        """
        if isinstance(uvw, numpy.ndarray):
            UVW = uvw
        else:
            UVW = numpy.array(uvw)

        T = self.to_cartesians_transform(as_array=True)
        XYZ = UVW @ T

        if as_array:
            return XYZ
        else:
            return XYZ.tolist()

    def to_cartesians_transform(self, as_array=False):
        """Matrix to convert fractional coordinates to Cartesian.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        as_array : bool = False
            Whether to return a numpy array or Python lists

        Returns
        -------
        transform : [N][float*3] or ndarray
            The transformation matrix
        """
        a, b, c, alpha, beta, gamma = self.parameters

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sg = sin(gamma)

        V = a * b * c * math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)
        # Transpose of ...
        # [a, b * cg, c * cb],
        # [0, b * sg, c * (ca - cb * cg) / sg],
        # [0, 0, V / (a * b * sg)]
        T = [
                [
                    a,
                    0,
                    0
                ],
                [
                    b * cg,
                    b * sg,
                    0
                ],
                [
                    c * cb,
                    c * (ca - cb * cg) / sg,
                    V / (a * b * sg)
                ]
            ]  # yapf: disable

        if as_array:
            return numpy.array(T)
        else:
            return T

    def to_fractionals(self, xyz, as_array=False):
        """Convert Cartesian coordinates to fractional.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        xyz : [N][3*float] or ndarray
            The Cartesian coordinates.

        Returns
        -------
        uvw : [N][float*3] or ndarray
            The ractional coordinates.
        """
        if isinstance(xyz, numpy.ndarray):
            XYZ = xyz
        else:
            XYZ = numpy.array(xyz)

        T = self.to_fractionals_transform(as_array=True)
        UVW = XYZ @ T

        if as_array:
            return UVW
        else:
            return UVW.tolist()

    def to_fractionals_transform(self, as_array=False):
        """Matrix to convert Cartesian coordinates to fractional.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        as_array : bool = False
            Whether to return a numpy array or Python lists

        Returns
        -------
        transform : [N][float*3] or ndarray
            The transformation matrix
        """
        a, b, c, alpha, beta, gamma = self.parameters

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sg = sin(gamma)

        V = a * b * c * math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)
        # Transpose...
        # [1 / a, -cg / (a * sg), b * c * (ca * cg - cb) / (V * sg)],
        # [0, 1 / (b * sg), a * c * (cb * cg - ca) / (V * sg)],
        # [0, 0, a * b * sg / V]
        T = [
                [
                    1 / a,
                    0,
                    0
                ],
                [
                    -cg / (a * sg),
                    1 / (b * sg),
                    0
                ],
                [
                    b * c * (ca * cg - cb) / (V * sg),
                    a * c * (cb * cg - ca) / (V * sg),
                    a * b * sg / V
                ]
            ]  # yapf: disable

        if as_array:
            return numpy.array(T)
        else:
            return T

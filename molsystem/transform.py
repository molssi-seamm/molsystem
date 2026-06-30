# -*- coding: utf-8 -*-

"""Rigid-body geometric transforms (translation, rotation) of a configuration.

The transforms operate on Cartesian coordinates in Å, modify the configuration
in place, and return ``self`` so that calls can be chained, e.g.::

    configuration.rotate(R, about="com").translate([0.0, 0.0, 5.0])

They are supported only for molecular (non-periodic) systems; calling them on a
periodic configuration raises a ``ValueError``. (Rigid-body moves of a periodic
system would also have to transform the unit cell to be meaningful, which is a
separate problem.)

The module also provides helpers for constructing rotation matrices
(:func:`rotation_matrix_from_axis_angle`, :func:`rotation_matrix_from_quaternion`,
:func:`random_rotation_matrix`).
"""

import math

import numpy as np


def rotation_matrix_from_axis_angle(axis, angle, degrees=True):
    """Build a 3x3 rotation matrix from an axis and angle (Rodrigues' formula).

    Parameters
    ----------
    axis : [float]*3
        The axis to rotate about. Need not be normalized.
    angle : float
        The angle to rotate through.
    degrees : bool = True
        Whether ``angle`` is in degrees (the default) or radians.

    Returns
    -------
    numpy.ndarray
        A 3x3 rotation matrix.
    """
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm == 0.0:
        raise ValueError("The rotation axis must be a non-zero vector.")
    x, y, z = axis / norm

    if degrees:
        angle = math.radians(angle)
    c = math.cos(angle)
    s = math.sin(angle)
    C = 1.0 - c

    return np.array(
        [
            [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
        ]
    )


def rotation_matrix_from_quaternion(quaternion):
    """Build a 3x3 rotation matrix from a quaternion (w, x, y, z).

    Parameters
    ----------
    quaternion : [float]*4
        The quaternion as (w, x, y, z). Need not be normalized.

    Returns
    -------
    numpy.ndarray
        A 3x3 rotation matrix.
    """
    q = np.asarray(quaternion, dtype=float)
    norm = np.linalg.norm(q)
    if norm == 0.0:
        raise ValueError("The quaternion must be a non-zero vector.")
    w, x, y, z = q / norm

    return np.array(
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - z * w),
                2.0 * (x * z + y * w),
            ],
            [
                2.0 * (x * y + z * w),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - x * w),
            ],
            [
                2.0 * (x * z - y * w),
                2.0 * (y * z + x * w),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    )


def random_rotation_matrix(rng=None):
    """A rotation matrix sampled uniformly from SO(3).

    Uses Shoemake's algorithm to draw a uniformly distributed unit quaternion.

    Parameters
    ----------
    rng : numpy.random.Generator = None
        The random-number generator to use. If None, a fresh default generator
        is created. Pass a seeded generator for reproducible orientations.

    Returns
    -------
    numpy.ndarray
        A 3x3 rotation matrix.
    """
    if rng is None:
        rng = np.random.default_rng()
    u1, u2, u3 = rng.random(3)
    r1 = math.sqrt(1.0 - u1)
    r2 = math.sqrt(u1)
    two_pi = 2.0 * math.pi
    quaternion = (
        r2 * math.cos(two_pi * u3),  # w
        r1 * math.sin(two_pi * u2),  # x
        r1 * math.cos(two_pi * u2),  # y
        r2 * math.sin(two_pi * u3),  # z
    )
    return rotation_matrix_from_quaternion(quaternion)


class TransformMixin:
    """A mixin giving a configuration rigid-body transform operations."""

    def center_of_mass(self, atoms=None, mass_weighted=True):
        """The center of mass (or geometric centroid) of the configuration.

        Parameters
        ----------
        atoms : [int] = None
            Optional 0-based positional indices selecting a subset of atoms
            (e.g. one fragment of a complex). By default all atoms are used.
        mass_weighted : bool = True
            If True (the default) return the mass-weighted center of mass; if
            False return the unweighted geometric centroid.

        Returns
        -------
        numpy.ndarray
            The center of mass as a length-3 array, in Å.
        """
        xyz = self.atoms.get_coordinates(fractionals=False, as_array=True)
        xyz = np.asarray(xyz, dtype=float)
        if atoms is not None:
            xyz = xyz[list(atoms)]

        if mass_weighted:
            masses = np.asarray(self.atoms.atomic_masses, dtype=float)
            if atoms is not None:
                masses = masses[list(atoms)]
            total = masses.sum()
            if total == 0.0:
                raise ValueError("The total mass is zero; cannot form a COM.")
            return (masses[:, np.newaxis] * xyz).sum(axis=0) / total
        else:
            return xyz.mean(axis=0)

    def centroid(self, atoms=None):
        """The (unweighted) geometric centroid of the configuration.

        A convenience wrapper for ``center_of_mass(mass_weighted=False)``.
        """
        return self.center_of_mass(atoms=atoms, mass_weighted=False)

    def translate(self, displacement, atoms=None):
        """Translate the configuration (or a subset of atoms) in place.

        Parameters
        ----------
        displacement : [float]*3
            The displacement vector, in Å.
        atoms : [int] = None
            Optional 0-based positional indices selecting a subset of atoms.
            By default all atoms are moved.

        Returns
        -------
        self
        """
        self._require_nonperiodic("translate")
        d = np.asarray(displacement, dtype=float).reshape(3)
        return self._apply_to_coordinates(lambda c: c + d, atoms=atoms)

    def rotate(self, rotation, about="com", atoms=None):
        """Rotate the configuration (or a subset of atoms) in place.

        Parameters
        ----------
        rotation : array-like
            A 3x3 rotation matrix, or a 4x4 matrix whose upper-left 3x3 block is
            used. The helpers :func:`rotation_matrix_from_axis_angle`,
            :func:`rotation_matrix_from_quaternion`, and
            :func:`random_rotation_matrix` build suitable matrices.
        about : str or [float]*3 = "com"
            The point to rotate about: "com" (center of mass), "centroid"
            (geometric center), "origin", or an explicit length-3 point in Å.
            When a subset of atoms is given, "com"/"centroid" refer to that
            subset.
        atoms : [int] = None
            Optional 0-based positional indices selecting a subset of atoms.
            By default the whole configuration is rotated.

        Returns
        -------
        self
        """
        self._require_nonperiodic("rotate")
        R = np.asarray(rotation, dtype=float)
        if R.shape == (4, 4):
            R = R[:3, :3]
        if R.shape != (3, 3):
            raise ValueError(
                f"The rotation must be a 3x3 or 4x4 matrix, not shape {R.shape}."
            )

        center = self._resolve_point(about, atoms=atoms)
        return self._apply_to_coordinates(
            lambda c: (c - center) @ R.T + center, atoms=atoms
        )

    def transform(self, matrix, atoms=None):
        """Apply a general affine transform in place (RDKit ``TransformMol`` form).

        The transform of a point ``p`` is ``R @ p + t``, where ``R`` and ``t``
        are taken from a 4x4 homogeneous matrix (or ``R`` from a 3x3 matrix with
        ``t`` = 0). This matches the convention used by RDKit's
        ``AllChem.TransformMol``, so a matrix obtained from RDKit alignment can
        be applied directly.

        Parameters
        ----------
        matrix : array-like
            A 3x3 (rotation only) or 4x4 (rotation + translation) matrix.
        atoms : [int] = None
            Optional 0-based positional indices selecting a subset of atoms.

        Returns
        -------
        self
        """
        self._require_nonperiodic("transform")
        M = np.asarray(matrix, dtype=float)
        if M.shape == (3, 3):
            R = M
            t = np.zeros(3)
        elif M.shape == (4, 4):
            R = M[:3, :3]
            t = M[:3, 3]
        else:
            raise ValueError(
                f"The transform must be a 3x3 or 4x4 matrix, not shape {M.shape}."
            )
        return self._apply_to_coordinates(lambda c: c @ R.T + t, atoms=atoms)

    def _require_nonperiodic(self, operation):
        """Raise if this configuration is periodic."""
        if self.periodicity != 0:
            raise ValueError(
                f"'{operation}' is only supported for non-periodic (molecular) "
                f"systems; this configuration has periodicity {self.periodicity}."
            )

    def _resolve_point(self, point, atoms=None):
        """Resolve a 'com'/'centroid'/'origin'/explicit point to an array."""
        if isinstance(point, str):
            if point == "com":
                return self.center_of_mass(atoms=atoms, mass_weighted=True)
            elif point == "centroid":
                return self.center_of_mass(atoms=atoms, mass_weighted=False)
            elif point == "origin":
                return np.zeros(3)
            else:
                raise ValueError(
                    f"Unknown reference point '{point}'; expected 'com', "
                    "'centroid', 'origin', or a length-3 point."
                )
        return np.asarray(point, dtype=float).reshape(3)

    def _apply_to_coordinates(self, func, atoms=None):
        """Apply ``func`` to the Cartesian coordinates and store the result.

        ``func`` takes and returns an (n, 3) array of Cartesian coordinates in Å.
        If ``atoms`` is given, only those rows are passed through ``func`` and
        updated; the rest are left unchanged.
        """
        xyz = self.atoms.get_coordinates(fractionals=False, as_array=True)
        xyz = np.asarray(xyz, dtype=float)
        if atoms is None:
            xyz = func(xyz)
        else:
            idx = list(atoms)
            xyz[idx] = func(xyz[idx])
        self.atoms.set_coordinates(xyz, fractionals=False)
        return self

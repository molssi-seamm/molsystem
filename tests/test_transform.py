#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the rigid-body transforms on a configuration."""

import math

import numpy as np
import pytest

from molsystem import (
    random_rotation_matrix,
    rotation_matrix_from_axis_angle,
    rotation_matrix_from_quaternion,
)


def _distance_matrix(xyz):
    """Pairwise distance matrix, for checking rigid-body invariance."""
    diff = xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :]
    return np.linalg.norm(diff, axis=-1)


# --------------------------------------------------------------------------- #
# Rotation-matrix helpers
# --------------------------------------------------------------------------- #


def test_axis_angle_90_about_z():
    """90 degrees about z maps (1, 0, 0) -> (0, 1, 0)."""
    R = rotation_matrix_from_axis_angle([0, 0, 1], 90)
    rotated = R @ np.array([1.0, 0.0, 0.0])
    assert np.allclose(rotated, [0.0, 1.0, 0.0], atol=1.0e-12)


def test_axis_angle_radians():
    """Degrees=False is honored."""
    R = rotation_matrix_from_axis_angle([0, 0, 1], math.pi / 2, degrees=False)
    rotated = R @ np.array([1.0, 0.0, 0.0])
    assert np.allclose(rotated, [0.0, 1.0, 0.0], atol=1.0e-12)


def test_axis_angle_unnormalized_axis():
    """The axis need not be normalized."""
    R1 = rotation_matrix_from_axis_angle([0, 0, 1], 37)
    R2 = rotation_matrix_from_axis_angle([0, 0, 5], 37)
    assert np.allclose(R1, R2)


def test_axis_angle_zero_axis_raises():
    with pytest.raises(ValueError):
        rotation_matrix_from_axis_angle([0, 0, 0], 90)


def test_quaternion_identity():
    """The identity quaternion gives the identity matrix."""
    R = rotation_matrix_from_quaternion([1, 0, 0, 0])
    assert np.allclose(R, np.eye(3))


def test_quaternion_is_proper_rotation():
    R = rotation_matrix_from_quaternion([0.5, 0.5, 0.5, 0.5])
    assert np.allclose(R @ R.T, np.eye(3), atol=1.0e-12)
    assert math.isclose(np.linalg.det(R), 1.0, abs_tol=1.0e-12)


def test_random_rotation_is_proper_and_reproducible():
    rng = np.random.default_rng(42)
    R = random_rotation_matrix(rng)
    assert np.allclose(R @ R.T, np.eye(3), atol=1.0e-12)
    assert math.isclose(np.linalg.det(R), 1.0, abs_tol=1.0e-12)
    # Same seed -> same matrix
    assert np.allclose(R, random_rotation_matrix(np.random.default_rng(42)))


# --------------------------------------------------------------------------- #
# center_of_mass / centroid
# --------------------------------------------------------------------------- #


def test_center_of_mass_vs_centroid(H2O):
    com = H2O.center_of_mass()
    centroid = H2O.centroid()
    # COM is pulled toward the heavy oxygen at the origin, so it sits below
    # the geometric centroid in z.
    assert com[2] < centroid[2]
    # By symmetry both lie on the z axis.
    assert np.allclose(com[:2], [0.0, 0.0], atol=1.0e-12)
    assert np.allclose(centroid[:2], [0.0, 0.0], atol=1.0e-12)


def test_centroid_is_mean(H2O):
    xyz = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    assert np.allclose(H2O.centroid(), xyz.mean(axis=0))


def test_center_of_mass_subset(AceticAcid):
    """COM of a subset equals the COM computed by hand for those atoms."""
    xyz = np.asarray(AceticAcid.atoms.get_coordinates(fractionals=False, as_array=True))
    masses = np.asarray(AceticAcid.atoms.atomic_masses)
    idx = [0, 4, 5]  # the two carbons and an oxygen
    expected = (masses[idx, np.newaxis] * xyz[idx]).sum(axis=0) / masses[idx].sum()
    assert np.allclose(AceticAcid.center_of_mass(atoms=idx), expected)


# --------------------------------------------------------------------------- #
# translate
# --------------------------------------------------------------------------- #


def test_translate_shifts_com(AceticAcid):
    com0 = AceticAcid.center_of_mass()
    d = np.array([1.0, -2.0, 3.5])
    AceticAcid.translate(d)
    assert np.allclose(AceticAcid.center_of_mass(), com0 + d)


def test_translate_preserves_internal_geometry(AceticAcid):
    d0 = _distance_matrix(
        np.asarray(AceticAcid.atoms.get_coordinates(fractionals=False, as_array=True))
    )
    AceticAcid.translate([3.0, 0.0, 0.0])
    d1 = _distance_matrix(
        np.asarray(AceticAcid.atoms.get_coordinates(fractionals=False, as_array=True))
    )
    assert np.allclose(d0, d1)


def test_translate_returns_self(AceticAcid):
    assert AceticAcid.translate([0.0, 0.0, 1.0]) is AceticAcid


def test_translate_subset(H2O):
    """Translating only the oxygen moves just that atom."""
    xyz0 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    H2O.translate([0.0, 0.0, 10.0], atoms=[0])
    xyz1 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    assert np.allclose(xyz1[0], xyz0[0] + [0.0, 0.0, 10.0])
    assert np.allclose(xyz1[1:], xyz0[1:])


# --------------------------------------------------------------------------- #
# rotate
# --------------------------------------------------------------------------- #


def test_rotate_about_com_preserves_com_and_shape(AceticAcid):
    com0 = AceticAcid.center_of_mass()
    d0 = _distance_matrix(
        np.asarray(AceticAcid.atoms.get_coordinates(fractionals=False, as_array=True))
    )
    R = random_rotation_matrix(np.random.default_rng(7))
    AceticAcid.rotate(R, about="com")
    assert np.allclose(AceticAcid.center_of_mass(), com0)
    d1 = _distance_matrix(
        np.asarray(AceticAcid.atoms.get_coordinates(fractionals=False, as_array=True))
    )
    assert np.allclose(d0, d1)


def test_rotate_about_origin_known(H2O):
    """90 degrees about z, about the origin: (x, 0, z) -> (0, x, z)."""
    xyz0 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    R = rotation_matrix_from_axis_angle([0, 0, 1], 90)
    H2O.rotate(R, about="origin")
    xyz1 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    expected = xyz0 @ R.T
    assert np.allclose(xyz1, expected)
    # Oxygen sits at the origin, so it is unmoved.
    assert np.allclose(xyz1[0], [0.0, 0.0, 0.0], atol=1.0e-12)


def test_rotate_accepts_4x4(H2O):
    """A 4x4 matrix uses only its rotation block in rotate()."""
    R = rotation_matrix_from_axis_angle([1, 1, 0], 33)
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = [9.0, 9.0, 9.0]  # translation block must be ignored by rotate()
    a = H2O.center_of_mass()
    H2O.rotate(M, about="com")
    # rotation about the COM leaves the COM fixed; translation block ignored
    assert np.allclose(H2O.center_of_mass(), a)


def test_rotate_bad_shape_raises(H2O):
    with pytest.raises(ValueError):
        H2O.rotate(np.ones((2, 2)))


def test_rotate_bad_about_raises(H2O):
    with pytest.raises(ValueError):
        H2O.rotate(np.eye(3), about="nonsense")


# --------------------------------------------------------------------------- #
# transform (RDKit TransformMol convention: p' = R @ p + t)
# --------------------------------------------------------------------------- #


def test_transform_matches_rotate_then_translate(AceticAcid, H2O):
    """A 4x4 transform equals rotating about the origin then translating."""
    R = rotation_matrix_from_axis_angle([0, 1, 0], 50)
    t = np.array([2.0, -1.0, 0.5])

    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = t

    # H2O and AceticAcid fixtures are independent configurations; build the
    # reference by hand from H2O's coordinates.
    xyz0 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    expected = xyz0 @ R.T + t

    H2O.transform(M)
    xyz1 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    assert np.allclose(xyz1, expected)


def test_transform_3x3_is_pure_rotation(H2O):
    xyz0 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    R = rotation_matrix_from_axis_angle([0, 0, 1], 90)
    H2O.transform(R)
    xyz1 = np.asarray(H2O.atoms.get_coordinates(fractionals=False, as_array=True))
    assert np.allclose(xyz1, xyz0 @ R.T)


def test_transform_bad_shape_raises(H2O):
    with pytest.raises(ValueError):
        H2O.transform(np.ones((2, 3)))


# --------------------------------------------------------------------------- #
# Periodic systems are disallowed for all transforms
# --------------------------------------------------------------------------- #


def test_translate_rejects_periodic(copper):
    with pytest.raises(ValueError):
        copper.translate([1.0, 0.0, 0.0])


def test_rotate_rejects_periodic(copper):
    with pytest.raises(ValueError):
        copper.rotate(np.eye(3))


def test_transform_rejects_periodic(copper):
    with pytest.raises(ValueError):
        copper.transform(np.eye(3))

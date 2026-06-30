#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for System.create_combined_configuration()."""

import math

import numpy as np
import pytest

from molsystem import rotation_matrix_from_axis_angle


def _add_water(configuration):
    """Populate a configuration with a single TIP3P-geometry water."""
    r0 = 0.9572
    theta0 = 104.52
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, z, z]
    atno = [8, 1, 1]

    ids = configuration.atoms.append(x=X, y=Y, z=Z, atno=atno)
    configuration.bonds.append(i=[ids[0], ids[0]], j=[ids[1], ids[2]], bondorder=[1, 1])
    return configuration


@pytest.fixture()
def two_waters(empty_db):
    """Two water monomers (in separate systems) plus an empty destination, one db."""
    db = empty_db
    cA = db.create_system(name="A").create_configuration(name="A")
    _add_water(cA)
    cB = db.create_system(name="B").create_configuration(name="B")
    _add_water(cB)
    dimers = db.create_system(name="dimers")
    return dimers, cA, cB


def _n_bonds(configuration):
    return len(configuration.bonds.get_column_data("i"))


def test_combine_counts(two_waters):
    dimers, cA, cB = two_waters
    new = dimers.create_combined_configuration([cA, cB])
    assert new.n_atoms == cA.n_atoms + cB.n_atoms == 6
    assert _n_bonds(new) == _n_bonds(cA) + _n_bonds(cB) == 4


def test_combine_atomic_numbers_concatenate(two_waters):
    dimers, cA, cB = two_waters
    new = dimers.create_combined_configuration([cA, cB])
    assert list(new.atoms.atomic_numbers) == (
        list(cA.atoms.atomic_numbers) + list(cB.atoms.atomic_numbers)
    )


def test_combine_places_with_transform(two_waters):
    dimers, cA, cB = two_waters
    shift = np.array([0.0, 0.0, 5.0])
    M = np.eye(4)
    M[:3, 3] = shift

    new = dimers.create_combined_configuration([cA, cB], transforms=[None, M])

    xyz = np.asarray(new.atoms.get_coordinates(fractionals=False, as_array=True))
    xyzA = np.asarray(cA.atoms.get_coordinates(fractionals=False, as_array=True))
    xyzB = np.asarray(cB.atoms.get_coordinates(fractionals=False, as_array=True))

    assert np.allclose(xyz[:3], xyzA)
    assert np.allclose(xyz[3:], xyzB + shift)


def test_combine_with_rotation_preserves_monomer_geometry(two_waters):
    dimers, cA, cB = two_waters
    R = rotation_matrix_from_axis_angle([1.0, 2.0, 3.0], 57.0)

    new = dimers.create_combined_configuration([cA, cB], transforms=[None, R])

    xyz = np.asarray(new.atoms.get_coordinates(fractionals=False, as_array=True))
    xyzB = np.asarray(cB.atoms.get_coordinates(fractionals=False, as_array=True))

    # Internal O-H distances of monomer B are unchanged by the rotation.
    def oh(coords):
        return np.linalg.norm(coords[1:] - coords[0], axis=1)

    assert np.allclose(sorted(oh(xyz[3:])), sorted(oh(xyzB)))


def test_combine_bonds_stay_within_monomers(two_waters):
    dimers, cA, cB = two_waters
    new = dimers.create_combined_configuration([cA, cB])

    ids = new.atoms.ids
    position = {atom_id: k for k, atom_id in enumerate(ids)}
    nA = cA.n_atoms

    i_ids = new.bonds.get_column_data("i")
    j_ids = new.bonds.get_column_data("j")
    assert len(i_ids) == 4
    for a, b in zip(i_ids, j_ids):
        # Both ends of every bond belong to the same monomer block.
        assert (position[a] < nA) == (position[b] < nA)


def test_combine_does_not_modify_sources(two_waters):
    dimers, cA, cB = two_waters
    before = np.asarray(cB.atoms.get_coordinates(fractionals=False, as_array=True))
    M = np.eye(4)
    M[:3, 3] = [10.0, 0.0, 0.0]

    dimers.create_combined_configuration([cA, cB], transforms=[None, M])

    after = np.asarray(cB.atoms.get_coordinates(fractionals=False, as_array=True))
    assert np.allclose(before, after)


def test_combine_sums_charge(empty_db):
    db = empty_db
    cA = db.create_system(name="A").create_configuration(name="A")
    _add_water(cA)
    cA.charge = -1
    cB = db.create_system(name="B").create_configuration(name="B")
    _add_water(cB)
    cB.charge = 2
    dest = db.create_system(name="dest")

    new = dest.create_combined_configuration([cA, cB])
    assert new.charge == 1


def test_combine_three_configurations(two_waters):
    """A trimer: combining three monomers works (n-body assembly)."""
    dimers, cA, cB = two_waters
    cC = cA.system_db.create_system(name="C").create_configuration(name="C")
    _add_water(cC)

    new = dimers.create_combined_configuration([cA, cB, cC])
    assert new.n_atoms == 9
    assert _n_bonds(new) == 6


def test_combine_rejects_periodic(two_waters):
    dimers, cA, cB = two_waters
    periodic = cA.system_db.create_system(name="P").create_configuration(name="P")
    periodic.periodicity = 3
    periodic.coordinate_system = "fractional"
    periodic.cell.parameters = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    periodic.atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol=["Ar"])

    with pytest.raises(ValueError):
        dimers.create_combined_configuration([cA, periodic])


def test_combine_transform_length_mismatch_raises(two_waters):
    dimers, cA, cB = two_waters
    with pytest.raises(ValueError):
        dimers.create_combined_configuration([cA, cB], transforms=[None])


def test_combine_empty_raises(two_waters):
    dimers, cA, cB = two_waters
    with pytest.raises(ValueError):
        dimers.create_combined_configuration([])

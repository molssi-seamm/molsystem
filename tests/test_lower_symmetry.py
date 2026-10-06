#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Configuration.lower_symmetry: in place and from another configuration."""

import pytest  # noqa: F401


def pair(system, name="pair", coordinate_system="Fractional", filler=True):
    """A cubic cell (a = 5 Å) holding a bonded N2 pair, 1 Å apart along a.

    With ``filler``, another configuration's atoms come first, so that the
    pair's atom ids do not start at the beginning.
    """
    if filler:
        other = system.create_configuration(
            name="filler", periodicity=3, coordinate_system="Fractional"
        )
        other.cell.parameters = [5.0, 5.0, 5.0, 90.0, 90.0, 90.0]
        other.atoms.append(x=[0.1] * 3, y=[0.0] * 3, z=[0.0] * 3, symbol=["He"] * 3)
    configuration = system.create_configuration(
        name=name, periodicity=3, coordinate_system=coordinate_system
    )
    configuration.cell.parameters = [5.0, 5.0, 5.0, 90.0, 90.0, 90.0]
    if coordinate_system == "Cartesian":
        x = [1.0, 2.0]
    else:
        x = [0.2, 0.4]
    ids = configuration.atoms.append(x=x, y=[0.0, 0.0], z=[0.0, 0.0], symbol=["N"] * 2)
    configuration.bonds.append(i=[ids[0]], j=[ids[1]], bondorder=[3])
    return configuration


def test_in_place_with_bonds(system):
    configuration = pair(system)
    configuration.lower_symmetry()
    assert configuration.atoms.n_atoms == 2
    assert configuration.bonds.n_bonds == 1
    assert configuration.bonds.get_column_data("bondorder") == [3]
    assert configuration.bonds.get_lengths() == pytest.approx([1.0])


@pytest.mark.parametrize("coordinate_system", ["Fractional", "Cartesian"])
def test_from_another_configuration(system, coordinate_system):
    source = pair(system, coordinate_system=coordinate_system)
    before = source.atoms.get_coordinates(fractionals=False, as_array=True)
    copy = system.create_configuration(
        name="copy", periodicity=3, coordinate_system=coordinate_system
    )
    copy.lower_symmetry(other=source)
    # The copy has the atoms, where they were, and the bond
    assert copy.atoms.get_coordinates(
        fractionals=False, as_array=True
    ) == pytest.approx(before)
    assert copy.bonds.n_bonds == 1
    assert copy.bonds.get_lengths() == pytest.approx([1.0])
    # in atoms and bonds of its own; the source is untouched
    assert copy.atomset != source.atomset and copy.bondset != source.bondset
    assert source.atoms.n_atoms == 2 and source.bonds.n_bonds == 1
    assert source.atoms.get_coordinates(
        fractionals=False, as_array=True
    ) == pytest.approx(before)


def test_space_group_expanded(system):
    """Face-centred cubic with one asymmetric atom: four atoms in P1."""
    configuration = system.create_configuration(
        name="fcc", periodicity=3, coordinate_system="Fractional"
    )
    configuration.symmetry.group = "F m -3 m"
    configuration.cell.parameters = [4.0, 4.0, 4.0, 90.0, 90.0, 90.0]
    configuration.atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol=["Cu"])
    configuration.lower_symmetry()
    assert configuration.symmetry.group == "P 1"
    assert sorted(
        tuple(round(v, 6) for v in xyz)
        for xyz in configuration.atoms.get_coordinates(fractionals=True)
    ) == [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]


def test_bond_across_a_symmetry_operation(system):
    """P -1: one asymmetric N, bonded to its inverse; a pair and its bond in P1."""
    configuration = system.create_configuration(
        name="pm1", periodicity=3, coordinate_system="Fractional"
    )
    configuration.symmetry.group = "P -1"
    configuration.cell.parameters = [5.0, 5.0, 5.0, 90.0, 90.0, 90.0]
    ids = configuration.atoms.append(x=[0.1], y=[0.0], z=[0.0], symbol=["N"])
    configuration.bonds.append(
        i=[ids[0]], j=[ids[0]], bondorder=[3], symop1=["."], symop2=["2"]
    )
    configuration.lower_symmetry()
    assert configuration.symmetry.group == "P 1"
    assert configuration.atoms.n_atoms == 2
    assert configuration.bonds.n_bonds == 1
    assert configuration.bonds.get_column_data("bondorder") == [3]
    assert configuration.bonds.get_lengths() == pytest.approx([1.0])

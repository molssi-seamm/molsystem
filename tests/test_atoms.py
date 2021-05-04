#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pprint  # noqa: F401
import pytest  # noqa: F401

import molsystem  # noqa: F401

"""Tests for `molsystem` package."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]


def test_construction(configuration):
    """Simplest test that we can make an Atoms object"""
    atoms = configuration.atoms
    assert str(type(atoms)) == "<class 'molsystem.atoms._Atoms'>"


def test_keys(atoms):
    """Test the default keys in an Atoms object"""
    result = ["atno", "configuration", "id", "x", "y", "z"]

    assert sorted([*atoms.keys()]) == result


def test_n_atoms_empty(atoms):
    """Test how many atoms an empty object has"""
    assert atoms.n_atoms == 0


def test_append_one(atoms):
    """Test adding one atom"""
    with atoms as tmp:
        tmp.append(x=1.0, y=2.0, z=3.0, atno=[6])
    assert atoms.n_atoms == 1


def test_append_several(atoms):
    """Test adding several atoms"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert atoms.n_atoms == 3


def test_append_several_scalar(atoms):
    """Test adding several atoms with some scalar values"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=0.0, atno=6)
    assert atoms["x"] == x
    assert atoms["y"] == y
    assert atoms["z"] == [0.0, 0.0, 0.0]
    assert atoms["atno"] == [6, 6, 6]


def test_append_using_default(atoms):
    """Test adding several atoms"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z)

    assert atoms.n_atoms == 3
    assert atoms["atno"] == [None, None, None]


def test_append_error(atoms):
    """Test adding atoms with an error"""
    with pytest.raises(KeyError) as e:
        with atoms as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno, bad=99)
    assert atoms.n_atoms == 0
    assert atoms.configuration.version == 0
    assert str(e.value) == "'\"bad\" is not an attribute of the atoms.'"


def test_add_attribute(atoms):
    """Test adding an attribute"""
    result = ["atno", "configuration", "id", "name", "x", "y", "z"]
    with atoms as tmp:
        tmp.add_attribute("name")
    assert sorted([*atoms.keys()]) == result


def test_add_duplicate_attribute(atoms):
    """Test duplicate adding an attribute"""
    result = ["atno", "configuration", "id", "name", "x", "y", "z"]
    with atoms as tmp:
        tmp.add_attribute("name")
    with pytest.raises(RuntimeError) as e:
        with atoms as tmp:
            ver = tmp.configuration.version
            tmp.add_attribute("name")
    assert sorted([*atoms.keys()]) == result
    assert ver == 1
    assert atoms.configuration.version == 1
    assert str(e.value) == "_Table attribute 'name' is already defined!"


def test_add_coordinates_attribute(atoms):
    """Test adding an attribute"""
    result = ["atno", "configuration", "id", "spin", "x", "y", "z"]
    with atoms as tmp:
        tmp.add_attribute("spin", configuration_dependent=True)
    assert sorted([*atoms.keys()]) == result
    del atoms["spin"]
    result.remove("spin")
    assert sorted([*atoms.keys()]) == result


def test_add_attribute_with_values(atoms):
    """Test adding several atoms"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("name", values=["H1", "O", "H2"])
    assert np.array_equal(atoms["name"], ["H1", "O", "H2"])


def test_add_attribute_with_one_value(atoms):
    """Test adding several atoms"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("name", values=["H1"])
    assert np.array_equal(atoms["name"], ["H1", "H1", "H1"])


def test_add_attribute_with_one_value_not_list(atoms):
    """Test adding an attribute using a scalar value"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("spin", coltype="int", default=0, values=1)
    assert atoms["spin"] == [1, 1, 1]


def test_add_attribute_with_wrong_number_of_values(atoms):
    """Test adding an attribute using a the wrong number of values"""
    with pytest.raises(IndexError) as e:
        with atoms as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno)
            tmp.add_attribute("spin", coltype="int", default=0, values=[1, 2])
    assert atoms.n_atoms == 0
    assert str(e.value) == (
        "The number of values given, "
        '2, must be either 1, or the number of rows in "main"."atom": 3'
    )


def test_get_attribute_by_index(atoms):
    """Get a single value of an attribute by index"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert atoms["atno"][1] == 1


def test_get_attribute_by_slice(atoms):
    """Get several values using a slice"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert atoms["atno"][::2] == [8, 1]


def test_getting_x(atoms):
    """Test getting the x column of coordinates"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert atoms["x"] == [1.0, 2.0, 3.0]


def test_contains(atoms):
    """Test the __contains__ or 'in' functionalty"""
    assert "x" in atoms


def test_not_in(atoms):
    """Test the __contains__ or 'in' functionalty"""
    assert "abc" not in atoms


def test_deleting_column(atoms):
    """Test deleting a column"""
    with atoms as tmp:
        del tmp["atno"]
    assert sorted([*atoms.keys()]) == ["configuration", "id", "x", "y", "z"]


def test_set_column(atoms):
    """Test setting a column using a scalar"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with atoms as tmp:
        tmp["atno"] = 10

    assert atoms.n_atoms == 3
    assert atoms.configuration.version == 2
    assert atoms["atno"] == [10, 10, 10]


def test_set_column_with_array(atoms):
    """Test setting a column using an array"""
    values = [10.0, 11.0, 12.0]

    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with atoms as tmp:
        tmp["x"] = values

    assert atoms.n_atoms == 3
    assert atoms.configuration.version == 2
    assert atoms["x"] == values


def test_append_error_no_coordinates(atoms):
    """Test adding atoms without coordinates, raising an error"""
    with atoms as tmp:
        tmp.append(y=y, z=z, atno=atno)

    assert atoms.n_atoms == 3
    assert atoms["x"] == [None, None, None]


def test_append_error_invalid_column(atoms):
    """Test adding atoms without coordinates, raising an error"""
    with pytest.raises(KeyError) as e:
        with atoms as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno, junk=99)
    assert atoms.n_atoms == 0
    assert str(e.value) == "'\"junk\" is not an attribute of the atoms.'"


def test_append_error_invalid_length(atoms):
    """Test adding atoms without coordinates, raising an error"""
    with pytest.raises(IndexError) as e:
        with atoms as tmp:
            tmp.append(x=x, y=y, z=z, atno=[3, 4])
    assert atoms.n_atoms == 0
    assert str(e.value) == (
        'key "atno" has the wrong number of values, '
        "2. Should be 1 or the number of atoms (3)."
    )


def test_add_attribute_with_no_default(atoms):
    """Test adding an attribute with no default, then several atoms"""
    with atoms as tmp:
        tmp.add_attribute("new", coltype="float")
        tmp.append(x=x, y=y, z=z, atno=atno, new=[-1, -2, -3])

    assert atoms.n_atoms == 3
    assert atoms["new"] == [-1, -2, -3]


def test_equality(two_configurations):
    """Test whether two Atoms objects are equal"""
    configuration1, configuration2 = two_configurations
    atoms1 = configuration1.atoms
    atoms2 = configuration2.atoms

    with atoms1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with atoms2 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    assert atoms1.n_atoms == 3
    assert atoms1 == atoms2


def test_inequality(two_configurations):
    """Test whether two Atoms objects are equal"""
    configuration1, configuration2 = two_configurations
    atoms1 = configuration1.atoms
    atoms2 = configuration2.atoms

    with atoms1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with atoms2 as tmp:
        atno2 = list(atno)
        atno2[2] = 13
        tmp.append(x=x, y=y, z=z, atno=atno2)

    assert atoms1.n_atoms == 3
    assert atoms1 != atoms2


def test_str(atoms):
    """Test string representation of atoms object"""
    string_rep = """\
   atno  configuration    x    y    z
1     8              1  1.0  4.0  7.0
2     1              1  2.0  5.0  8.0
3     1              1  3.0  6.0  9.0"""  # noqa: W291

    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    if str(atoms) != string_rep:
        print(str(atoms))

    assert str(atoms) == string_rep


def test_repr(atoms):
    """Test representation of atoms object"""
    string_rep = """\
   atno  configuration    x    y    z
1     8              1  1.0  4.0  7.0
2     1              1  2.0  5.0  8.0
3     1              1  3.0  6.0  9.0"""  # noqa: W291

    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    if repr(atoms) != string_rep:
        print(repr(atoms))

    assert repr(atoms) == string_rep


def test_equals(atoms):
    """Test shallow copy"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    atoms2 = atoms

    assert atoms == atoms2


def test_equals_with_change(atoms):
    """Test (shallow) copy and changing copy"""
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    atoms2 = atoms
    atoms2["atno"][0] = 22

    assert atoms == atoms2


def test_coordinates(AceticAcid):
    """Test getting coordinates."""
    configuration = AceticAcid

    xs = [1.0797, 0.5782, 0.7209, 0.7052, 0.5713, -0.1323, 0.9757, 2.1724]
    ys = [0.0181, 3.1376, -0.6736, -0.3143, 1.3899, 1.7142, 2.2970, 0.0161]
    zs = [-0.0184, 0.2813, -0.7859, 0.9529, -0.3161, -1.2568, 0.5919, -0.0306]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    xyz = configuration.atoms.coordinates
    assert np.allclose(xyz, xyz0)


def test_set_coordinates(AceticAcid):
    """Test setting coordinates."""
    configuration = AceticAcid

    xs = [1.08, 0.58, 0.72, 0.71, 0.57, -0.13, 0.98, 2.17]
    ys = [0.02, 3.14, -0.67, -0.31, 1.39, 1.71, 2.30, 0.02]
    zs = [-0.02, 0.28, -0.79, 0.95, -0.32, -1.26, 0.59, -0.03]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    configuration.atoms.set_coordinates(xyz0)

    xyz = configuration.atoms.coordinates
    assert np.allclose(xyz, xyz0)


def test_selected_atoms(AceticAcid):
    """Test getting the number of selected atoms."""
    configuration = AceticAcid

    assert configuration.atoms.get_n_atoms("atno", "==", 6) == 2
    x = []
    for row in configuration.atoms.atoms("atno", "==", 6):
        x.append(row["x"])
    assert x == [1.0797, 0.5713]


def test_periodic_coordinates(vanadium):
    """Test getting coordinates."""
    configuration = vanadium

    xs = [0.0, 0.5]
    ys = [0.0, 0.5]
    zs = [0.0, 0.5]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    xyz = configuration.atoms.coordinates
    assert np.allclose(xyz, xyz0)


def test_periodic_coordinates_cartesians(vanadium):
    """Test getting coordinates in Cartesian coordinates."""
    configuration = vanadium

    xs = [0.0, 1.515]
    ys = [0.0, 1.515]
    zs = [0.0, 1.515]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    xyz = configuration.atoms.get_coordinates(fractionals=False)
    assert np.allclose(xyz, xyz0)


def test_set_periodic_coordinates(vanadium):
    """Test setting coordinates."""
    configuration = vanadium

    xs = [0.1, 0.6]
    ys = [0.1, 0.6]
    zs = [0.1, 0.6]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    configuration.atoms.set_coordinates(xyz0)

    xyz = configuration.atoms.coordinates
    assert np.allclose(xyz, xyz0)


def test_set_periodic_coordinates_cartesians(vanadium):
    """Test setting coordinates in Cartesian coordinates."""
    configuration = vanadium

    xs = [0.303, 1.818]
    ys = [0.303, 1.818]
    zs = [0.303, 1.818]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    configuration.atoms.set_coordinates(xyz0, fractionals=False)

    xyz = configuration.atoms.get_coordinates(fractionals=False)
    assert np.allclose(xyz, xyz0)

    xs = [0.1, 0.6]
    ys = [0.1, 0.6]
    zs = [0.1, 0.6]
    xyz0 = [[x, y, z] for x, y, z in zip(xs, ys, zs)]

    xyz = configuration.atoms.coordinates
    assert np.allclose(xyz, xyz0)


def test_remove_atoms(copper):
    """Test removing one atom from FCC copper"""
    configuration = copper

    assert configuration.atoms.n_atoms == 4

    ids = configuration.atoms.ids
    first = [ids[0]]

    with configuration.atoms as tmp:
        tmp.delete(atoms=first)

    assert configuration.atoms.n_atoms == 3

    xyz = configuration.atoms.coordinates

    assert xyz == [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]

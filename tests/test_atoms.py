#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy

import numpy as np
import pytest  # noqa: F401

import molsystem  # noqa: F401
"""Tests for `molsystem` package."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]

xa = np.array(x)
ya = np.array(y)
za = np.array(z)
atnoa = np.array(atno)


def test_construction():
    """Simplest test that we can make an Atoms object"""
    system = None
    atoms = molsystem.Atoms(system)
    assert str(type(atoms)) == "<class 'molsystem.atoms.Atoms'>"


def test_keys():
    """Test the default keys in an Atoms object"""
    system = None
    atoms = molsystem.Atoms(system)
    assert sorted([*atoms.keys()]) == ['atno', 'x', 'y', 'z']


def test_n_atoms_empty():
    """Test how many atoms an empty object has"""
    system = None
    atoms = molsystem.Atoms(system)
    assert atoms.n_atoms == 0


def test_append_one():
    """Test adding one atom"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=1.0, y=2.0, z=3.0, atno=[6])
    assert atoms.n_atoms == 1


def test_append_several():
    """Test adding several atoms"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert atoms.n_atoms == 3


def test_append_several_scalar():
    """Test adding several atoms with some scalar values"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=x, y=y, z=0.0, atno=6)
    assert (
        np.array_equal(atoms['x'], x) and np.array_equal(atoms['y'], y) and
        np.array_equal(atoms['z'], [0.0, 0.0, 0.0]) and
        np.array_equal(atoms['atno'], [6, 6, 6])
    )


def test_append_ndarrays():
    """Test adding several atoms us ndarrays for data"""
    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert atoms.n_atoms == 3


def test_append_ndarrays_using_default():
    """Test adding several atoms"""
    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za)

    assert (atoms.n_atoms == 3 and np.array_equal(atoms['atno'], [-1, -1, -1]))


def test_append_ndarrays_3x100():
    """Test adding 3 atoms 100 times, to force allocation"""
    system = None
    atoms = molsystem.Atoms(system)
    xaa = np.array(xa)
    yaa = np.array(ya)
    zaa = np.array(za)
    atnoaa = np.array(atnoa)
    with atoms as tmp:
        for i in range(0, 100):
            tmp.append(x=xaa, y=yaa, z=zaa, atno=atnoaa)
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert atoms.n_atoms == 300 and atoms.version == 1


def test_append_ndarrays_3x10_one_by_one():
    """Test adding 3 atoms 10 times separately, checking changes"""
    system = None
    atoms = molsystem.Atoms(system)
    xaa = np.array(xa)
    yaa = np.array(ya)
    zaa = np.array(za)
    atnoaa = np.array(atnoa)
    for i in range(0, 10):
        with atoms as tmp:
            tmp.append(x=xaa, y=yaa, z=zaa, atno=atnoaa)
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert atoms.n_atoms == 30 and atoms.version == 10


def test_trim():
    """Test adding several atoms and then reclaiming empty space"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        free_space = tmp.free
        tmp.trim()

    assert (
        len(atoms) == 3 and atoms.free == 0 and
        free_space == (molsystem.Atoms.allocate_min - 3)
    )


def test_append_error():
    """Test adding atoms with an error"""
    system = None
    atoms = molsystem.Atoms(system)
    try:
        with atoms as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
            raise RuntimeError()
    except RuntimeError:
        pass
    assert atoms.n_atoms == 0 and atoms.version == 0


def test_add_attribute():
    """Test adding an attribute"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.add_attribute('name')
    assert sorted([*atoms.keys()]) == ['atno', 'name', 'x', 'y', 'z']


def test_add_duplicate_attribute():
    """Test duplicate adding an attribute"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.add_attribute('name')
    try:
        with atoms as tmp:
            ver = tmp.version
            tmp.add_attribute('name')
    except RuntimeError as e:
        err = str(e)
    assert (
        sorted([*atoms.keys()]) == ['atno', 'name', 'x', 'y', 'z'] and
        ver == 1 and atoms.version == 1 and
        err == "Atoms attribute 'name' is already defined!"
    )


def test_define_duplicate_attribute():
    """Test duplicate definition of an attribute"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.add_attribute('name')
    try:
        with atoms as tmp:
            tmp.define_attribute('name', coltype=np.object, default='')
    except RuntimeError as e:
        err = str(e)
    assert (err == "Atoms attribute 'name' is already defined!")


def test_add_attribute_with_different_type():
    """Add an attribute giving a different type"""
    system = None
    atoms = molsystem.Atoms(system)
    try:
        with atoms as tmp:
            tmp.add_attribute('name', coltype=np.float64)
    except ValueError as e:
        err = str(e)
    assert (
        err == (
            "Column type should be '<class 'object'>', "
            "not '<class 'numpy.float64'>'"
        )
    )


def test_add_attribute_with_different_default():
    """Add an attribute giving a different default"""
    system = None
    atoms = molsystem.Atoms(system)
    try:
        with atoms as tmp:
            tmp.add_attribute('name', default='unknown')
    except ValueError as e:
        err = str(e)
    assert (err == "Default should be '', not 'unknown'")


def test_add_attribute_with_values():
    """Test adding several atoms"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('name', values=['H1', 'O', 'H2'])
    assert np.array_equal(atoms['name'], ['H1', 'O', 'H2'])


def test_add_attribute_with_one_value():
    """Test adding several atoms"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('name', values=['H1'])
    assert np.array_equal(atoms['name'], ['H1', 'H1', 'H1'])


def test_add_attribute_with_one_value_not_list():
    """Test adding an attribute using a scalar value"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('spin', coltype=np.int8, default=0, values=1)
    assert np.array_equal(atoms['spin'], [1, 1, 1])


def test_add_attribute_with_wrong_number_of_values():
    """Test adding an attribute using a the wrong number of values"""
    system = None
    atoms = molsystem.Atoms(system)
    try:
        with atoms as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
            tmp.add_attribute(
                'spin', coltype=np.int8, default=0, values=[1, 2]
            )
    except IndexError as e:
        err = str(e)

    assert (
        atoms.n_atoms == 0 and err == (
            "The number of values given, "
            "2, must be either 1, or the number of atoms: 3"
        )
    )


def test_get_attribute_by_index():
    """Get a single value of an attribute by index"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert atoms['atno'][1] == 1


def test_get_attribute_by_slice():
    """Get several values using a slice"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert np.array_equal(atoms['atno'][::2], [8, 1])


def test_getting_x():
    """Test getting the x column of coordinates"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert np.array_equal(atoms['x'], [1.0, 2.0, 3.0])


def test_contains():
    """Test the __contains__ or 'in' functionalty"""
    system = None
    atoms = molsystem.Atoms(system)
    assert 'x' in atoms


def test_not_in():
    """Test the __contains__ or 'in' functionalty"""
    system = None
    atoms = molsystem.Atoms(system)
    assert 'abc' not in atoms


def test_deleting_column():
    """Test deleting a column"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        del tmp['atno']
    assert sorted([*atoms.keys()]) == ['x', 'y', 'z']


def test_set_column():
    """Test setting a column using a scalar"""
    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with atoms as tmp:
        tmp['atno'] = 10

    assert (
        atoms.n_atoms == 3 and atoms.version == 2 and
        np.array_equal(atoms['atno'], [10, 10, 10])
    )


def test_set_column_with_array():
    """Test setting a column using an array"""
    system = None
    atoms = molsystem.Atoms(system)

    values = [10.0, 11.0, 12.0]

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with atoms as tmp:
        tmp['x'] = values

    assert (
        atoms.n_atoms == 3 and atoms.version == 2 and
        np.array_equal(atoms['x'], values)
    )


def test_append_error_no_coordinates():
    """Test adding atoms without coordinates, raising an error"""
    system = None
    atoms = molsystem.Atoms(system)

    try:
        with atoms as tmp:
            tmp.append(y=ya, z=za, atno=atnoa)
    except KeyError as e:
        err = str(e)

    assert (atoms.n_atoms == 0 and err == "'The coordinates are required!'")


def test_append_error_invalid_column():
    """Test adding atoms without coordinates, raising an error"""
    system = None
    atoms = molsystem.Atoms(system)

    try:
        with atoms as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa, junk=3)
    except KeyError as e:
        err = str(e)

    assert (
        atoms.n_atoms == 0 and
        err == '\'"junk" is not an attribute of the atoms!\''
    )


def test_append_error_invalid_length():
    """Test adding atoms without coordinates, raising an error"""
    system = None
    atoms = molsystem.Atoms(system)

    try:
        with atoms as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=[3, 4])
    except IndexError as e:
        err = str(e)

    assert (
        atoms.n_atoms == 0 and err == (
            'key "atno" has the wrong number of values, '
            '2. Should be 1 or the number of atoms (3).'
        )
    )


def test_add_attribute_with_no_default():
    """Test adding an attribute with no default, then several atoms"""
    system = None
    atoms = molsystem.Atoms(system)
    with atoms as tmp:
        tmp.add_attribute('new', coltype=np.float64)
        tmp.append(x=xa, y=ya, z=za, atno=atnoa, new=[-1, -2, -3])

    assert (atoms.n_atoms == 3 and np.array_equal(atoms['new'], [-1, -2, -3]))


def test_add_attribute_with_no_default_error():
    """Test adding an attribute with no default, then several atoms"""
    system = None
    atoms = molsystem.Atoms(system)
    try:
        with atoms as tmp:
            tmp.add_attribute('new', coltype=np.float64)
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    except KeyError as e:
        err = str(e)
    assert (
        atoms.n_atoms == 0 and err == (
            '"There is no default for attribute '
            "'new'. You must supply a value\""
        )
    )


def test_set_free():
    """Test using atoms.free = xxxx to set free space"""
    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.free = 1000

    assert (atoms.n_atoms == 3 and atoms.free == 1000 and len(atoms) == 1003)


def test_equality():
    """Test whether two Atoms objects are equal"""
    system = None
    atoms1 = molsystem.Atoms(system)
    atoms2 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with atoms2 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert atoms1.n_atoms == 3 and atoms1 == atoms2


def test_inequality():
    """Test whether two Atoms objects are equal"""
    system = None
    atoms1 = molsystem.Atoms(system)
    atoms2 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with atoms2 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp['atno'][2] = 13

    assert atoms1.n_atoms == 3 and atoms1 != atoms2


def test_str():
    """Test string representation of atoms object"""
    string_rep = """\
       x    y    z  atno
uid                     
0    1.0  4.0  7.0   8.0
1    2.0  5.0  8.0   1.0
2    3.0  6.0  9.0   1.0"""  # noqa: W291

    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert str(atoms) == string_rep


def test_repr():
    """Test representation of atoms object"""
    string_rep = """\
       x    y    z  atno
uid                     
0    1.0  4.0  7.0   8.0
1    2.0  5.0  8.0   1.0
2    3.0  6.0  9.0   1.0"""  # noqa: W291

    system = None
    atoms = molsystem.Atoms(system)

    with atoms as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    print(repr(atoms))

    assert repr(atoms) == string_rep


def test_copy_constructor():
    """Test copy constructor"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = molsystem.Atoms(atoms1)

    assert atoms1 == atoms2


def test_copy_constructor_with_change():
    """Test copy constructor and changing copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = molsystem.Atoms(atoms1)
    atoms2['atno'][0] = 22

    assert atoms1 != atoms2


def test_equals():
    """Test shallow copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = atoms1

    assert atoms1 == atoms2


def test_equals_with_change():
    """Test (shallow) copy and changing copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = atoms1
    atoms2['atno'][0] = 22

    assert atoms1 == atoms2


def test_copy():
    """Test copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = copy.copy(atoms1)

    assert atoms1 == atoms2


def test_copy_with_change():
    """Test copy and changing copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = copy.copy(atoms1)
    atoms2['atno'][0] = 22

    assert atoms1 == atoms2


def test_deep_copy():
    """Test copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = copy.deepcopy(atoms1)

    assert atoms1 == atoms2


def test_deep_copy_with_change():
    """Test deepcopy and changing copy"""
    system = None
    atoms1 = molsystem.Atoms(system)

    with atoms1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    atoms2 = copy.deepcopy(atoms1)
    atoms2['atno'][0] = 22

    assert atoms1 != atoms2


def test_coordinate_type():
    atoms = molsystem.Atoms()

    with atoms as tmp:
        tmp.coordinate_type = 'f'
    assert atoms.coordinate_type == 'fractional'

    with atoms as tmp:
        tmp.coordinate_type = 'cart'
    assert atoms.coordinate_type == 'Cartesian'


def test_invalid_coordinate_value():
    atoms = molsystem.Atoms()

    with pytest.raises(
        ValueError,
        match=("The coordinate_type must be 'Cartesian' or 'fractional', .*")
    ):
        atoms.coordinate_type = 'bad!'


def test_invalid_coordinate_type():
    atoms = molsystem.Atoms()

    with pytest.raises(
        ValueError,
        match=("The coordinate_type must be 'Cartesian' or 'fractional', .*")
    ):
        atoms.coordinate_type = 1

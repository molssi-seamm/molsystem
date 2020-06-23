#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `bonds` in the `molsystem` package."""

import copy

import numpy as np
import pytest  # noqa: F401

from molsystem import Bonds

#         C  H  H  H  C =O  O  H
i_atom = [0, 0, 0, 0, 4, 4, 6]
j_atom = [1, 2, 3, 4, 5, 6, 7]
order = [1, 1, 1, 1, 1, 2, 1]


@pytest.fixture()
def bonds():
    """A bonds object initialize with acetic acid
    """
    bonds = Bonds()
    bonds.append(i=i_atom, j=j_atom, order=order)
    return bonds


@pytest.fixture(scope='session')
def bonds_session():
    """A bonds object initialized with acetic acid
    """
    bonds = Bonds()
    bonds.append(i=i_atom, j=j_atom, order=order)
    return bonds


def test_construction():
    """Simplest test that we can make an Bonds object"""
    system = None
    bonds = Bonds(system)
    assert str(type(bonds)) == "<class 'molsystem.bonds.Bonds'>"


def test_keys():
    """Test the default keys in an Bonds object"""
    system = None
    bonds = Bonds(system)
    assert sorted([*bonds.keys()]) == ['i', 'j', 'order']


def test_n_bonds_empty():
    """Test how many bonds an empty object has"""
    system = None
    bonds = Bonds(system)
    assert bonds.n_bonds == 0


def test_append_one():
    """Test adding one bond"""
    system = None
    bonds = Bonds(system)
    with bonds as tmp:
        tmp.append(i=(0,), j=(1,))
    assert bonds.n_bonds == 1


def test_append_several(bonds):
    """Test adding several bonds"""
    assert bonds.n_bonds == 7


def test_append_scalar_i():
    """Test adding bonds with one being a scalar"""
    system = None
    bonds = Bonds(system)
    with bonds as tmp:
        tmp.append(i=0, j=j_atom[0:4])
    assert bonds.n_bonds == 4


def test_append_scalar_j():
    """Test adding bonds with one being a scalar"""
    system = None
    bonds = Bonds(system)
    with bonds as tmp:
        tmp.append(j=0, i=j_atom[0:4])
    assert bonds.n_bonds == 4


def test_append_ndarrays():
    """Test adding several bonds us ndarrays for data"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    j = np.array([1, 2, 3, 4])
    order = np.array([1, 1, 1, 1])

    with bonds as tmp:
        tmp.append(i=i, j=j, order=order)

    assert bonds.n_bonds == 4


def test_append_ndarrays_using_default():
    """Test adding several bonds"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    j = np.array([1, 2, 3, 4])

    with bonds as tmp:
        tmp.append(i=i, j=j)

    assert (
        bonds.n_bonds == 4 and np.array_equal(bonds['order'], [1, 1, 1, 1])
    )


def test_append_ndarrays_4x100():
    """Test adding 4 bonds 100 times, to force allocation"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    j = np.array([1, 2, 3, 4])
    order = np.array([1, 1, 1, 1])

    with bonds as tmp:
        for _ in range(0, 100):
            tmp.append(i=i, j=j, order=order)
            i += 5
            j += 5
    assert bonds.n_bonds == 400 and bonds.version == 1


def test_append_ndarrays_4x10_one_by_one():
    """Test adding 4 bonds 10 times separately, checking changes"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    j = np.array([1, 2, 3, 4])
    order = np.array([1, 1, 1, 1])

    for _ in range(0, 10):
        with bonds as tmp:
            tmp.append(i=i, j=j, order=order)
            i += 5
            j += 5

    assert bonds.n_bonds == 40 and bonds.version == 10


def test_trim():
    """Test adding several bonds and then reclaiming empty space"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    j = np.array([1, 2, 3, 4])
    order = np.array([1, 1, 1, 1])

    with bonds as tmp:
        tmp.append(i=i, j=j, order=order)
        free_space = tmp.free
        tmp.trim()

    assert (
        len(bonds) == 4 and bonds.free == 0 and
        free_space == (Bonds.allocate_min - 4)
    )


def test_append_error():
    """Test adding bonds with an error"""
    system = None
    bonds = Bonds(system)

    i = np.array([0, 0, 0, 0])
    order = np.array([1, 1, 1, 1])

    try:
        with bonds as tmp:
            tmp.append(i=i, order=order)
    except KeyError:
        pass
    assert bonds.n_bonds == 0 and bonds.version == 0


def test_add_attribute(bonds):
    """Test adding an attribute"""
    with bonds as tmp:
        tmp.add_attribute('name')
    assert sorted([*bonds.keys()]) == ['i', 'j', 'name', 'order']


def test_add_duplicate_attribute(bonds):
    """Test duplicate adding an attribute"""
    with bonds as tmp:
        tmp.add_attribute('name')
    try:
        with bonds as tmp:
            ver = tmp.version
            tmp.add_attribute('name')
    except RuntimeError as e:
        err = str(e)
    assert (
        sorted([*bonds.keys()]) == ['i', 'j', 'name', 'order'] and ver == 1 and
        bonds.version == 1 and
        err == "Bonds attribute 'name' is already defined!"
    )


def test_define_duplicate_attribute(bonds):
    """Test duplicate definition of an attribute"""
    with bonds as tmp:
        tmp.add_attribute('name')
    try:
        with bonds as tmp:
            tmp.define_attribute('name', coltype=np.object, default='')
    except RuntimeError as e:
        err = str(e)
    assert (err == "Bonds attribute 'name' is already defined!")


def test_add_attribute_with_values(bonds):
    """Test adding several bonds"""
    names = ['CH1', 'CH2', 'CH3', 'CC', 'C=O', 'CO', 'OH']
    with bonds as tmp:
        tmp.add_attribute('name', values=names)
    assert np.array_equal(bonds['name'], names)


def test_add_attribute_with_one_value(bonds):
    """Test adding several bonds"""
    with bonds as tmp:
        tmp.add_attribute('name', values=['b'])
    assert np.array_equal(bonds['name'], 7 * ['b'])


def test_add_attribute_with_one_value_not_list(bonds):
    """Test adding an attribute using a scalar value"""
    with bonds as tmp:
        tmp.add_attribute('spin', coltype=np.int8, default=0, values=1)
    assert np.array_equal(bonds['spin'], 7 * [1])


def test_add_attribute_with_wrong_number_of_values(bonds):
    """Test adding an attribute using a the wrong number of values"""
    try:
        with bonds as tmp:
            tmp.add_attribute(
                'spin', coltype=np.int8, default=0, values=[1, 2]
            )
    except IndexError as e:
        err = str(e)

    assert (
        bonds.n_bonds == 7 and err == (
            "The number of values given, "
            "2, must be either 1, or the number of bonds: 7"
        )
    )


def test_get_attribute_by_index(bonds_session):
    """Get a single value of an attribute by index"""
    assert bonds_session['j'][1] == 2


def test_get_attribute_by_slice(bonds_session):
    """Get several values using a slice"""
    assert np.array_equal(bonds_session['i'][::2], [0, 0, 4, 6])


def test_contains(bonds_session):
    """Test the __contains__ or 'in' functionalty"""
    assert 'order' in bonds_session


def test_not_in(bonds_session):
    """Test the __contains__ or 'in' functionalty"""
    assert 'abc' not in bonds_session


def test_deleting_column(bonds):
    """Test deleting a column"""
    with bonds as tmp:
        del tmp['order']
    assert sorted([*bonds.keys()]) == ['i', 'j']


def test_set_column(bonds):
    """Test setting a column using a scalar"""
    with bonds as tmp:
        tmp['order'] = 2

    assert (
        bonds.n_bonds == 7 and bonds.version == 1 and
        np.array_equal(bonds['order'], 7 * [2])
    )


def test_set_column_with_array(bonds):
    """Test setting a column using an array"""
    values = [1, 2, 3, 1, 2, 3, 1]

    with bonds as tmp:
        tmp['order'] = values

    assert (
        bonds.n_bonds == 7 and bonds.version == 1 and
        np.array_equal(bonds['order'], values)
    )


def test_append_error_no_atom(bonds):
    """Test adding bonds without an atom, raising an error"""
    try:
        with bonds as tmp:
            tmp.append(i=(2,))
    except KeyError as e:
        err = str(e)

    assert (bonds.n_bonds == 7 and err == "'The atoms i & j are required!'")


def test_append_error_invalid_length():
    """Test adding bonds with different numbers of items"""
    system = None
    bonds = Bonds(system)

    try:
        with bonds as tmp:
            tmp.append(i=[0, 0, 0], j=[1, 2])
    except IndexError as e:
        err = str(e)
        print(err)

    assert (
        bonds.n_bonds == 0 and err == (
            'key "j" has the wrong number of values, 2. Should be 1 or the '
            'number of values in i, 3.'
        )
    )


def test_add_attribute_with_no_default():
    """Test adding an attribute with no default, then several bonds"""
    system = None
    bonds = Bonds(system)
    with bonds as tmp:
        tmp.add_attribute('new', coltype=np.float64)
        tmp.append(i=[0, 0, 0, 0], j=[1, 2, 3, 4], new=[-1, -2, -3, -4])

    assert (
        bonds.n_bonds == 4 and np.array_equal(bonds['new'], [-1, -2, -3, -4])
    )


def test_add_attribute_with_no_default_error():
    """Test adding an attribute with no default, then several bonds"""
    system = None
    bonds = Bonds(system)
    try:
        with bonds as tmp:
            tmp.add_attribute('new', coltype=np.float64)
            tmp.append(i=[0, 0, 0, 0], j=[1, 2, 3, 4])
    except KeyError as e:
        err = str(e)
    assert (
        bonds.n_bonds == 0 and err == (
            '"There is no default for attribute '
            "'new'. You must supply a value\""
        )
    )


def test_set_free(bonds):
    """Test using bonds.free = xxxx to set free space"""
    with bonds as tmp:
        tmp.free = 1000

    assert (bonds.n_bonds == 7 and bonds.free == 1000 and len(bonds) == 1007)


def test_equality(bonds, bonds_session):
    """Test whether two Bonds objects are equal"""
    assert bonds.n_bonds == 7 and bonds == bonds_session


def test_inequality(bonds, bonds_session):
    """Test whether two Bonds objects are equal"""
    with bonds as tmp:
        tmp['i'][2] = 13

    assert bonds.n_bonds == 7 and bonds != bonds_session


def test_str(bonds_session):
    """Test string representation of bonds object"""
    string_rep = """\
       i    j  order
uid                 
0    0.0  1.0    1.0
1    0.0  2.0    1.0
2    0.0  3.0    1.0
3    0.0  4.0    1.0
4    4.0  5.0    1.0
5    4.0  6.0    2.0
6    6.0  7.0    1.0"""  # noqa: W291

    # print(bonds_session)
    assert str(bonds_session) == string_rep


def test_repr(bonds_session):
    """Test representation of bonds object"""
    string_rep = """\
       i    j  order
uid                 
0    0.0  1.0    1.0
1    0.0  2.0    1.0
2    0.0  3.0    1.0
3    0.0  4.0    1.0
4    4.0  5.0    1.0
5    4.0  6.0    2.0
6    6.0  7.0    1.0"""  # noqa: W291

    # print(repr(bonds_session))
    assert repr(bonds_session) == string_rep


def test_copy_constructor(bonds_session):
    """Test copy constructor"""
    bonds = Bonds(bonds_session)

    assert bonds == bonds_session


def test_copy_constructor_with_change(bonds_session):
    """Test copy constructor and changing copy"""
    bonds = Bonds(bonds_session)

    with bonds as tmp:
        tmp['order'][0] = 2

    assert bonds != bonds_session


def test_equals(bonds_session):
    """Test shallow copy"""
    bonds = bonds_session

    assert bonds == bonds_session


def test_equals_with_change(bonds):
    """Test (shallow) copy and changing copy"""
    bonds2 = bonds
    bonds2['order'][0] = 2

    assert bonds == bonds2


def test_copy(bonds_session):
    """Test copy"""

    bonds = copy.copy(bonds_session)

    assert bonds == bonds_session


def test_copy_with_change(bonds):
    """Test copy and changing copy"""
    bonds2 = copy.copy(bonds)
    bonds2['order'][0] = 2

    assert bonds == bonds2


def test_deep_copy(bonds_session):
    """Test deep copy"""
    bonds2 = copy.deepcopy(bonds_session)

    assert bonds_session == bonds2


def test_deep_copy_with_change(bonds):
    """Test deepcopy and changing copy"""
    bonds2 = copy.deepcopy(bonds)
    bonds2['order'][0] = 2

    assert bonds != bonds2


def test_invalid_type(bonds):
    """Try using a bad type, e.g. a string."""
    with pytest.raises(
        TypeError, match=("'i' and 'j', the atoms indices, must be integers")
    ):
        with bonds as tmp:
            tmp.append(i=('i1', 'i2'), j=(3, 4))

    assert bonds.version == 0


def test_invalid_type2(bonds):
    """Try using a bad type, e.g. a string."""
    with pytest.raises(
        TypeError, match=("'i' and 'j', the atoms indices, must be integers")
    ):
        with bonds as tmp:
            tmp.append(i=(1, 2), j=('j1', 'j2'))

    assert bonds.version == 0

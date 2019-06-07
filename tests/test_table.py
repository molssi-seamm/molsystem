#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy

import numpy as np
import pytest  # noqa: F401

import molsystem  # noqa: F401
"""Tests for the Table class."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]

xa = np.array(x)
ya = np.array(y)
za = np.array(z)
atnoa = np.array(atno)

template = molsystem.Table()
template._private['attributes'] = {
    'x': {
        'type': np.float64,
        'default': 0.0,
    },
    'y': {
        'type': np.float64,
        'default': 0.0,
    },
    'z': {
        'type': np.float64,
        'default': 0.0,
    },
    'atno': {
        'type': np.int8,
        'default': -1,
    },
    'name': {
        'type': np.object,
        'default': '',
    },
}

for column in ('x', 'y', 'z', 'atno'):
    template.add_attribute(column)


def test_construction():
    """Simplest test that we can make an Table object"""
    table = molsystem.Table(template)
    assert str(type(table)) == "<class 'molsystem.table.Table'>"


def test_keys():
    """Test the default keys in an Table object"""
    table = molsystem.Table(template)
    assert sorted([*table.keys()]) == ['atno', 'x', 'y', 'z']


def test_n_rows_empty():
    """Test how many rows an empty object has"""
    table = molsystem.Table(template)
    assert table.n_rows == 0


def test_append_one():
    """Test adding one row"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=1.0, y=2.0, z=3.0, atno=[6])
    assert table.n_rows == 1


def test_append_several():
    """Test adding several rows"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert table.n_rows == 3


def test_append_several_scalar():
    """Test adding several table with some scalar values"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=x, y=y, z=0.0, atno=6)
    assert (
        np.array_equal(table['x'], x) and np.array_equal(table['y'], y) and
        np.array_equal(table['z'], [0.0, 0.0, 0.0]) and
        np.array_equal(table['atno'], [6, 6, 6])
    )


def test_append_several_scalar_first():
    """Test adding several table with a scalar value first"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(atno=6, x=x, y=y, z=0.0)
    assert (
        np.array_equal(table['x'], x) and np.array_equal(table['y'], y) and
        np.array_equal(table['z'], [0.0, 0.0, 0.0]) and
        np.array_equal(table['atno'], [6, 6, 6])
    )


def test_append_ndarrays():
    """Test adding several rows using ndarrays for data"""
    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert table.n_rows == 3


def test_append_ndarrays_using_default():
    """Test adding several rows using default values"""
    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za)

    assert (table.n_rows == 3 and np.array_equal(table['atno'], [-1, -1, -1]))


def test_append_ndarrays_3x100():
    """Test adding 3 rows 100 times, to force allocation"""
    table = molsystem.Table(template)
    xaa = np.array(xa)
    yaa = np.array(ya)
    zaa = np.array(za)
    atnoaa = np.array(atnoa)
    with table as tmp:
        for i in range(0, 100):
            tmp.append(x=xaa, y=yaa, z=zaa, atno=atnoaa)
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert table.n_rows == 300 and table.version == 1


def test_append_ndarrays_3x10_one_by_one():
    """Test adding 3 rows 10 times separately, checking changes"""
    table = molsystem.Table(template)
    xaa = np.array(xa)
    yaa = np.array(ya)
    zaa = np.array(za)
    atnoaa = np.array(atnoa)
    for i in range(0, 10):
        with table as tmp:
            tmp.append(x=xaa, y=yaa, z=zaa, atno=atnoaa)
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert table.n_rows == 30 and table.version == 10


def test_trim():
    """Test adding several rows and then reclaiming empty space"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        free_space = tmp.free
        tmp.trim()

    assert (
        len(table) == 3 and table.free == 0 and
        free_space == (molsystem.Table.allocate_min - 3)
    )


def test_append_error():
    """Test adding rows with an error"""
    table = molsystem.Table(template)
    try:
        with table as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
            raise RuntimeError()
    except RuntimeError:
        pass
    assert table.n_rows == 0 and table.version == 0


def test_add_attribute():
    """Test adding an attribute"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.add_attribute('name')
    assert sorted([*table.keys()]) == ['atno', 'name', 'x', 'y', 'z']


def test_add_duplicate_attribute():
    """Test duplicate adding an attribute"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.add_attribute('name')
    try:
        with table as tmp:
            ver = tmp.version
            tmp.add_attribute('name')
    except RuntimeError as e:
        err = str(e)
    assert (
        sorted([*table.keys()]) == ['atno', 'name', 'x', 'y', 'z'] and
        ver == 1 and table.version == 1 and
        err == "Table attribute 'name' is already defined!"
    )


def test_define_duplicate_attribute():
    """Test duplicate definition of an attribute"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.add_attribute('name')
    try:
        with table as tmp:
            tmp.define_attribute('name', coltype=np.object, default='')
    except RuntimeError as e:
        err = str(e)
    assert (err == "Table attribute 'name' is already defined!")


def test_add_attribute_with_different_type():
    """Add an attribute giving a different type"""
    table = molsystem.Table(template)
    try:
        with table as tmp:
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
    table = molsystem.Table(template)
    try:
        with table as tmp:
            tmp.add_attribute('name', default='unknown')
    except ValueError as e:
        err = str(e)
    assert (err == "Default should be '', not 'unknown'")


def test_add_attribute_with_values():
    """Test adding an attribute with values"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('name', values=['H1', 'O', 'H2'])
    assert np.array_equal(table['name'], ['H1', 'O', 'H2'])


def test_add_attribute_with_one_value():
    """Test adding an attribute with a list length 1"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('name', values=['H1'])
    assert np.array_equal(table['name'], ['H1', 'H1', 'H1'])


def test_add_attribute_with_one_value_not_list():
    """Test adding an attribute using a scalar value"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.add_attribute('spin', coltype=np.int8, default=0, values=1)
    assert np.array_equal(table['spin'], [1, 1, 1])


def test_add_attribute_with_wrong_number_of_values():
    """Test adding an attribute using the wrong number of values"""
    table = molsystem.Table(template)
    try:
        with table as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
            tmp.add_attribute(
                'spin', coltype=np.int8, default=0, values=[1, 2]
            )
    except IndexError as e:
        err = str(e)

    assert (
        table.n_rows == 0 and err == (
            "The number of values given, "
            "2, must be either 1, or the number of rows: 3"
        )
    )


def test_get_attribute_by_index():
    """Get a single value of an attribute by index"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert table['atno'][1] == 1


def test_get_attribute_by_slice():
    """Get several values using a slice"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert np.array_equal(table['atno'][::2], [8, 1])


def test_getting_x():
    """Test getting the x column of coordinates"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    assert np.array_equal(table['x'], [1.0, 2.0, 3.0])


def test_contains():
    """Test the __contains__ or 'in' functionalty"""
    table = molsystem.Table(template)
    assert 'x' in table


def test_not_in():
    """Test the __contains__ or 'in' functionalty"""
    table = molsystem.Table(template)
    assert 'abc' not in table


def test_deleting_column():
    """Test deleting a column"""
    table = molsystem.Table(template)
    with table as tmp:
        del tmp['atno']
    assert sorted([*table.keys()]) == ['x', 'y', 'z']


def test_no_change():
    """Test not making a change"""
    table = molsystem.Table(template)
    with table as tmp:
        ver = tmp.version
    assert ver == 0 and table.version == 0


def test_set_column():
    """Test setting a column using a scalar"""
    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with table as tmp:
        tmp['atno'] = 10

    assert (
        table.n_rows == 3 and table.version == 2 and
        np.array_equal(table['atno'], [10, 10, 10])
    )


def test_set_column_with_array():
    """Test setting a column using an array"""
    table = molsystem.Table(template)

    values = [10.0, 11.0, 12.0]

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with table as tmp:
        tmp['x'] = values

    assert (
        table.n_rows == 3 and table.version == 2 and
        np.array_equal(table['x'], values)
    )


def test_append_error_invalid_column():
    """Test adding rows with an unknown attribute, raising an error"""
    table = molsystem.Table(template)

    try:
        with table as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=atnoa, junk=3)
    except KeyError as e:
        err = str(e)

    assert (
        table.n_rows == 0 and
        err == '\'"junk" is not an attribute of the table!\''
    )


def test_append_error_invalid_length():
    """Test appending with the wrong number of values, raising an error"""
    table = molsystem.Table(template)

    try:
        with table as tmp:
            tmp.append(x=xa, y=ya, z=za, atno=[3, 4])
    except IndexError as e:
        err = str(e)

    assert (
        table.n_rows == 0 and err == (
            'key "atno" has the wrong number of values, '
            '2. Should be 1 or the number of rows (3).'
        )
    )


def test_add_attribute_with_no_default():
    """Test adding an attribute with no default, then several rows"""
    table = molsystem.Table(template)
    with table as tmp:
        tmp.add_attribute('new', coltype=np.float64)
        tmp.append(x=xa, y=ya, z=za, atno=atnoa, new=[-1, -2, -3])

    assert (table.n_rows == 3 and np.array_equal(table['new'], [-1, -2, -3]))


def test_add_attribute_with_no_default_error():
    """Test adding several rows with no value for an attribute with no default
    """
    table = molsystem.Table(template)
    try:
        with table as tmp:
            tmp.add_attribute('new', coltype=np.float64)
            tmp.append(x=xa, y=ya, z=za, atno=atnoa)
    except KeyError as e:
        err = str(e)
    assert (
        table.n_rows == 0 and err == (
            '"There is no default for attribute '
            "'new'. You must supply a value\""
        )
    )


def test_set_free():
    """Test using table.free = xxxx to set free space"""
    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp.free = 1000

    assert (table.n_rows == 3 and table.free == 1000 and len(table) == 1003)


def test_equality():
    """Test whether two Table objects are equal"""
    table1 = molsystem.Table(template)
    table2 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with table2 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert table1.n_rows == 3 and table1 == table2


def test_inequality():
    """Test whether two Table objects are equal"""
    table1 = molsystem.Table(template)
    table2 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    with table2 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)
        tmp['atno'][2] = 13

    assert table1.n_rows == 3 and table1 != table2


def test_str():
    """Test string representation of table object"""
    string_rep = """\
       x    y    z  atno
uid                     
0    1.0  4.0  7.0   8.0
1    2.0  5.0  8.0   1.0
2    3.0  6.0  9.0   1.0"""  # noqa: W291

    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    assert str(table) == string_rep


def test_repr():
    """Test representation of table object"""
    string_rep = """\
       x    y    z  atno
uid                     
0    1.0  4.0  7.0   8.0
1    2.0  5.0  8.0   1.0
2    3.0  6.0  9.0   1.0"""  # noqa: W291

    table = molsystem.Table(template)

    with table as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    print(repr(table))

    assert repr(table) == string_rep


def test_copy_constructor():
    """Test copy constructor"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = molsystem.Table(table1)

    assert table1 == table2


def test_copy_constructor_with_change():
    """Test copy constructor and changing copy"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = molsystem.Table(table1)
    table2['atno'][0] = 22

    assert table1 != table2


def test_equlas():
    """Test copy using equals"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = table1

    assert table1 == table2


def test_equals_with_change():
    """Test copy using equals and changing copy"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = table1
    table2['atno'][0] = 22

    assert table1 == table2


def test_copy_with_change():
    """Test copy and changing copy"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = copy.copy(table1)
    table2['atno'][0] = 22

    assert table1 == table2


def test_deep_copy():
    """Test copy"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = copy.deepcopy(table1)

    assert table1 == table2


def test_deep_copy_with_change():
    """Test deepcopy and changing copy"""
    table1 = molsystem.Table(template)

    with table1 as tmp:
        tmp.append(x=xa, y=ya, z=za, atno=atnoa)

    table2 = copy.deepcopy(table1)
    table2['atno'][0] = 22

    assert table1 != table2

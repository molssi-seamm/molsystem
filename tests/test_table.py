#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import pprint

import numpy as np
import pytest  # noqa: F401

import molsystem  # noqa: F401

"""Tests for the Table class."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]


def test_construction(simple_table):
    """Simplest test that we can make an Table object"""
    assert str(type(simple_table)) == "<class 'molsystem.table._Table'>"


def test_keys(simple_table):
    """Test the default keys in an Table object"""
    assert sorted([*simple_table.keys()]) == ["atno", "x", "y", "z"]


def test_n_rows_empty(simple_table):
    """Test how many rows an empty object has"""
    assert simple_table.n_rows == 0


def test_append_one(simple_table):
    """Test adding one row"""
    with simple_table as tmp:
        tmp.append(x=1.0, y=2.0, z=3.0, atno=[6])
    assert simple_table.n_rows == 1


def test_append_several(simple_table):
    """Test adding several rows"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert simple_table.n_rows == 3


def test_append_several_scalar(simple_table):
    """Test adding several simple_table with some scalar values"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=0.0, atno=6)
    assert (
        simple_table["x"] == x
        and simple_table["y"] == y
        and simple_table["z"] == [0.0, 0.0, 0.0]
        and simple_table["atno"] == [6, 6, 6]
    )


def test_append_several_scalar_first(simple_table):
    """Test adding several simple_table with a scalar value first"""
    with simple_table as tmp:
        tmp.append(atno=6, x=x, y=y, z=0.0)
    assert (
        simple_table["x"] == x
        and simple_table["y"] == y
        and simple_table["z"] == [0.0, 0.0, 0.0]
        and simple_table["atno"] == [6, 6, 6]
    )


def test_append_using_default(simple_table):
    """Test adding several rows using default values"""

    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z)

    assert simple_table.n_rows == 3 and simple_table["atno"] == [-1, -1, -1]


def test_append_ndarrays_3x100(simple_table):
    """Test adding 3 rows 100 times, to force allocation"""
    xaa = np.array(x)
    yaa = np.array(y)
    zaa = np.array(z)
    atnoaa = np.array(atno)
    with simple_table as tmp:
        for i in range(0, 100):
            tmp.append(
                x=xaa.tolist(), y=yaa.tolist(), z=zaa.tolist(), atno=atnoaa.tolist()
            )
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert simple_table.n_rows == 300


def test_append_ndarrays_3x10_one_by_one(simple_table):
    """Test adding 3 rows 10 times separately, checking changes"""
    xaa = np.array(x)
    yaa = np.array(y)
    zaa = np.array(z)
    atnoaa = np.array(atno)
    for i in range(0, 10):
        with simple_table as tmp:
            tmp.append(
                x=xaa.tolist(), y=yaa.tolist(), z=zaa.tolist(), atno=atnoaa.tolist()
            )
            xaa += 10.0
            yaa += 5.0
            zaa -= 3.0
    assert simple_table.n_rows == 30


def test_append_error(simple_table):
    """Test adding rows with an error"""
    try:
        with simple_table as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno)
            raise RuntimeError()
    except RuntimeError:
        pass
    assert simple_table.n_rows == 0


def test_add_attribute(simple_table):
    """Test adding an attribute"""
    with simple_table as tmp:
        tmp.add_attribute("name")
    assert sorted([*simple_table.keys()]) == ["atno", "name", "x", "y", "z"]


def test_add_duplicate_attribute(simple_table):
    """Test duplicate adding an attribute"""
    with simple_table as tmp:
        tmp.add_attribute("name")
    try:
        with simple_table as tmp:
            tmp.add_attribute("name")
    except RuntimeError as e:
        err = str(e)
    assert sorted([*simple_table.keys()]) == ["atno", "name", "x", "y", "z"]
    assert err == "_Table attribute 'name' is already defined!"


def test_add_attribute_with_values(simple_table):
    """Test adding an attribute with values"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("name", values=["H1", "O", "H2"])
    assert simple_table["name"] == ["H1", "O", "H2"]


def test_add_attribute_with_one_value(simple_table):
    """Test adding an attribute with a list length 1"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("name", values=["H1"])
    assert simple_table["name"] == ["H1", "H1", "H1"]


def test_add_attribute_with_one_value_not_list(simple_table):
    """Test adding an attribute using a scalar value"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("spin", coltype="int", default=0, values=1)
    assert simple_table["spin"] == [1, 1, 1]


def test_add_attribute_with_wrong_number_of_values(simple_table):
    """Test adding an attribute using the wrong number of values"""
    try:
        with simple_table as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno)
            tmp.add_attribute("spin", coltype="int", default=0, values=[1, 2])
    except IndexError as e:
        err = str(e)

    assert simple_table.n_rows == 0 and err == (
        "The number of values given, 2, must be either 1, or the number of"
        ' rows in "main"."table1": 3'
    )


def test_get_attribute_by_index(simple_table):
    """Get a single value of an attribute by index"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert simple_table["atno"][1] == 1


def test_get_attribute_by_slice(simple_table):
    """Get several values using a slice"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert simple_table["atno"][::2] == [8, 1]


def test_getting_x(simple_table):
    """Test getting the x column of coordinates"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
    assert simple_table["x"] == [1.0, 2.0, 3.0]


def test_contains(simple_table):
    """Test the __contains__ or 'in' functionalty"""
    assert "x" in simple_table


def test_not_in(simple_table):
    """Test the __contains__ or 'in' functionalty"""
    assert "abc" not in simple_table


def test_deleting_column(simple_table):
    """Test deleting a column"""
    with simple_table as tmp:
        del tmp["atno"]
    assert sorted([*simple_table.keys()]) == ["x", "y", "z"]


def test_set_column(simple_table):
    """Test setting a column using a scalar"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with simple_table as tmp:
        tmp["atno"] = 10

    assert simple_table.n_rows == 3
    assert simple_table["atno"] == [10, 10, 10]


def test_set_column_with_array(simple_table):
    """Test setting a column using an array"""
    values = [10.0, 11.0, 12.0]

    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with simple_table as tmp:
        tmp["x"] = values

    assert simple_table.n_rows == 3
    assert simple_table["x"] == values


def test_append_error_invalid_column(simple_table):
    """Test adding rows with an unknown attribute, raising an error"""
    try:
        with simple_table as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno, junk=3)
    except KeyError as e:
        err = str(e)

    assert (
        simple_table.n_rows == 0
        and err == '\'"junk" is not an attribute of the table "main"."table1"!\''
    )


def test_append_error_invalid_length(simple_table):
    """Test appending with the wrong number of values, raising an error"""
    try:
        with simple_table as tmp:
            tmp.append(x=x, y=y, z=z, atno=atno[0:2])
    except IndexError as e:
        err = str(e)

    assert simple_table.n_rows == 0 and err == (
        'key "atno" has the wrong number of values, '
        '2. Should be 1 or the number of rows in "main"."table1" (3).'
    )


def test_add_attribute_with_no_default(simple_table):
    """Test adding an attribute with no default, then several rows"""
    with simple_table as tmp:
        tmp.add_attribute("new", coltype="float")
        tmp.append(x=x, y=y, z=z, atno=atno, new=[-1, -2, -3])

    assert simple_table.n_rows == 3 and simple_table["new"] == [-1, -2, -3]


def test_add_attribute_with_no_default_error(simple_table):
    """Test adding a NOT NULL attribute with no default."""
    with pytest.raises(ValueError) as e:
        with simple_table as tmp:
            tmp.add_attribute("new", coltype="float", notnull=True)
    assert simple_table.n_rows == 0 and str(e.value) == (
        "Not null attributes must have defaults: new"
    )


def test_equality(two_tables):
    """Test whether two simple_table objects are equal"""
    table1, table2 = two_tables

    with table1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with table2 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    assert table1.n_rows == 3 and table1 == table2


def test_inequality(two_tables):
    """Test whether two simple_table objects are not equal"""
    table1, table2 = two_tables

    with table1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    with table2 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp["atno"][2] = 13

    assert table1.n_rows == 3 and table1 != table2


def test_str(simple_table):
    """Test string representation of simple_table object"""
    string_rep = """\
   atno    x    y    z
1     8  1.0  4.0  7.0
2     1  2.0  5.0  8.0
3     1  3.0  6.0  9.0"""  # noqa: W291

    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    assert str(simple_table) == string_rep


def test_repr(simple_table):
    """Test representation of simple_table object"""
    string_rep = """\
   atno    x    y    z
1     8  1.0  4.0  7.0
2     1  2.0  5.0  8.0
3     1  3.0  6.0  9.0"""  # noqa: W291

    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    assert repr(simple_table) == string_rep


def test_equals(simple_table):
    """Test copy using equals"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    simple_table2 = simple_table

    assert simple_table == simple_table2


def test_equals_with_change(simple_table):
    """Test copy using equals and changing copy"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    simple_table2 = simple_table

    atnos = simple_table2["atno"]
    atnos[2] = 22
    simple_table2["atno"] = atnos

    assert simple_table == simple_table2


def test_copy_with_change(simple_table):
    """Test copy and changing copy"""
    with simple_table as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    simple_table2 = copy.copy(simple_table)

    atnos = simple_table2["atno"]
    atnos[2] = 22
    simple_table2["atno"] = atnos

    assert simple_table == simple_table2


def test_diff(two_tables):
    """Test diffing two tables in the same database."""

    ref1 = {
        "changed": {1: {("atno", 8, 10), ("x", 1.0, 0.0)}},
        "columns in added rows": ["atno", "x", "y", "z"],
        "added": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {"columns changed": {"atno", "x"}, "rows added": 1},
    }

    ref2 = {
        "changed": {1: {("atno", 10, 8), ("x", 0.0, 1.0)}},
        "columns in deleted rows": ["atno", "x", "y", "z"],
        "deleted": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {"columns changed": {"atno", "x"}, "rows deleted": 1},
    }

    table1, table2 = two_tables

    with table1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    x1 = list(x)
    x1[0] = 0.0
    with table2 as tmp:
        tmp.append(x=x1, y=y, z=z, atno=[10, 1, 1])
        tmp.append(x=20.0, y=21.0, z=22.0, atno=12)

    diffs = table2.diff(table1)
    if diffs != ref1:
        pprint.pprint(diffs)
    assert diffs == ref1

    diffs = table1.diff(table2)
    if diffs != ref2:
        pprint.pprint(diffs)
    assert diffs == ref2


def test_diff_two_dbs(two_tables_in_two_dbs):
    """Test diffing two tables in the different databases."""

    ref1 = {
        "changed": {1: {("atno", 8, 10), ("x", 1.0, 0.0)}},
        "columns in added rows": ["atno", "x", "y", "z"],
        "added": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {"columns changed": {"atno", "x"}, "rows added": 1},
    }

    ref2 = {
        "changed": {1: {("atno", 10, 8), ("x", 0.0, 1.0)}},
        "columns in deleted rows": ["atno", "x", "y", "z"],
        "deleted": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {"columns changed": {"atno", "x"}, "rows deleted": 1},
    }

    table1, table2 = two_tables_in_two_dbs

    with table1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)

    x1 = list(x)
    x1[0] = 0.0
    with table2 as tmp:
        tmp.append(x=x1, y=y, z=z, atno=[10, 1, 1])
        tmp.append(x=20.0, y=21.0, z=22.0, atno=12)

    diffs = table2.diff(table1)
    if diffs != ref1:
        pprint.pprint(diffs)
    assert diffs == ref1

    diffs = table1.diff(table2)
    if diffs != ref2:
        pprint.pprint(diffs)
    assert diffs == ref2


def test_diff_add_columns(two_tables_in_two_dbs):
    """Test diffing two tables with different numbers of columns."""

    ref1 = {
        "columns deleted": {"spin"},
        "changed": {1: {("x", 1.0, 0.0), ("atno", 8, 10)}},
        "columns in added rows": ["atno", "x", "y", "z"],
        "added": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {
            "columns changed": {"atno", "x"},
            "columns deleted": {"spin"},
            "rows added": 1,
        },
    }  # yapf: disable

    ref2 = {
        "columns added": {"spin"},
        "changed": {1: {("x", 0.0, 1.0), ("atno", 10, 8)}},
        "columns in deleted rows": ["atno", "x", "y", "z"],
        "deleted": {4: (12, 20.0, 21.0, 22.0)},
        "summary": {
            "columns add": {"spin"},
            "columns changed": {"atno", "x"},
            "rows deleted": 1,
        },
    }  # yapf: disable

    table1, table2 = two_tables_in_two_dbs

    with table1 as tmp:
        tmp.append(x=x, y=y, z=z, atno=atno)
        tmp.add_attribute("spin", coltype="float", default=0.0)

    x1 = list(x)
    x1[0] = 0.0
    with table2 as tmp:
        tmp.append(x=x1, y=y, z=z, atno=[10, 1, 1])
        tmp.append(x=20.0, y=21.0, z=22.0, atno=12)

    diffs = table2.diff(table1)
    if diffs != ref1:
        pprint.pprint(diffs)
    assert diffs == ref1

    diffs = table1.diff(table2)
    if diffs != ref2:
        pprint.pprint(diffs)
    assert diffs == ref2

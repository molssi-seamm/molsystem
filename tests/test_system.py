#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the system class."""

import pytest  # noqa: F401


def test_construction(system):
    """Simplest test that we can make a System object"""
    assert str(type(system)) == "<class 'molsystem.system._System'>"


def test_version_empty(system):
    """Simplest test that we can make a System object"""
    assert system.version == 0


def test_create_table(system):
    """Test that we can create a table in the system."""
    with system as tmp:
        table = tmp.create_table('table1')
        # Create a column so the table exists
        table.add_attribute('data')
    assert system.version == 1
    assert 'table1' in system


def test_create_table_(system):
    """Test that we can create a table using [] syntax."""
    with system as tmp:
        table = tmp['table1']
        # Create a column so the table exists
        table.add_attribute('data')
    assert system.version == 1
    assert 'table1' in system


def test_delete_table(system_with_two_tables):
    """Test that we can delete a table."""
    system = system_with_two_tables
    with system as tmp:
        del tmp['table1']
    assert system.version == 1
    assert 'table2' in system
    assert 'table1' not in system


def test_elements(system):
    """Test that we can access the element table"""
    table = system['element']
    assert table.n_rows == 118

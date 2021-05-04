#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the system class."""

import pytest  # noqa: F401


def test_construction(system):
    """Simplest test that we can make a System object"""
    assert str(type(system)) == "<class 'molsystem.system._System'>"
    del system


def test_create_table(system):
    """Test that we can create a table in the system."""
    table = system.create_table("table1")
    # Create a column so the table exists
    table.add_attribute("data")
    assert "table1" in system


def test_create_table_element(system):
    """Test that we can create a table using [] syntax."""
    table = system["table1"]
    # Create a column so the table exists
    table.add_attribute("data")
    assert "table1" in system


def test_delete_table(system_with_two_tables):
    """Test that we can delete a table."""
    system = system_with_two_tables
    del system["table1"]
    assert "table2" in system
    assert "table1" not in system


def test_elements(system):
    """Test that we can access the element table"""
    table = system["element"]
    assert table.n_rows == 118


def test_create_new_system(AceticAcid):
    """Test creating a new, empty system."""
    system_db = AceticAcid.system_db
    assert system_db.n_systems == 1
    system_db.create_system(name="new system")
    assert system_db.n_systems == 2
    assert system_db.system_ids == [1, 2]


def test_create_new_system_vanadium(AceticAcid):
    """Test creating a new system, filling it for vanadium."""
    system_db = AceticAcid.system_db
    assert system_db.n_systems == 1

    system = system_db.create_system(name="BCC Vanadium")
    system.create_configuration(name="default")
    assert system_db.n_systems == 2
    assert system_db.system_ids == [1, 2]

    system_db.current_system_id = 2

    configuration = system_db.system.configuration

    configuration.periodicity = 3
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (3.03, 3.03, 3.03, 90, 90, 90)
    configuration.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol="V")

    formula, empirical_formula, Z = configuration.formula
    assert "".join(formula) == "V2"
    assert "".join(empirical_formula) == "V"
    assert Z == 2

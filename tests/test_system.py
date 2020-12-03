#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the system class."""

import pprint

import pytest  # noqa: F401


def test_construction(system):
    """Simplest test that we can make a System object"""
    assert str(type(system)) == "<class 'molsystem.system._System'>"
    del system


def test_version_empty(system):
    """Simplest test that we can make a System object"""
    assert system.version == 0
    del system


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


def test_formula(AceticAcid):
    """Test the formula generation."""
    formula, empirical_formula, Z = AceticAcid.formula()
    assert ''.join(formula) == 'C2H4O2'
    assert ''.join(empirical_formula) == 'CH2O'
    assert Z == 2


def test_formula_periodic(vanadium):
    """Test the formula generation."""
    formula, empirical_formula, Z = vanadium.formula()
    assert ''.join(formula) == 'V2'
    assert ''.join(empirical_formula) == 'V'
    assert Z == 2


def test_clear(CH3COOH_3H2O):
    """Test making subsets for the molecules."""
    result = [2, 3, 4, 5]

    system = CH3COOH_3H2O
    sids = system.create_molecule_subsets()
    if sids != result:
        pprint.pprint(sids)
    assert sids == result

    system.clear()

    assert system.n_atoms() == 0


def test_density(vanadium):
    """Test the density, and implicitly the mass and volume."""

    assert abs(vanadium.density() - 6.0817308915133) < 1.0e-06

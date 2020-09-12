#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cell` in the `molsystem` package."""

import pytest  # noqa: F401


def test_construction(system):
    """Simplest test that we can make a CellParameters object"""
    cell = system['cell']
    assert str(
        type(cell)
    ) == "<class 'molsystem.cell_parameters._CellParameters'>"


def test_periodic_system(system):
    """Test making a simple periodic system."""
    system.periodicity = 3
    system.coordinate_system = 'fractional'
    system.cell.set_cell(3.03, 3.03, 3.03, 90, 90, 90)
    system.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol='V')
    assert system.n_atoms() == 2 and system.version == 0


def test_periodic_context(system):
    """Test using context."""
    with system as tmp:
        tmp.periodicity = 3
        tmp.coordinate_system = 'fractional'
        tmp.cell.set_cell(3.03, 3.03, 3.03, 90, 90, 90)
        tmp.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol='V')
    assert system.n_atoms() == 2 and system.version == 1


def test_cell_object(vanadium):
    """Test getting the cell as a Cell object."""
    cell = vanadium.cell.cell()
    assert str(type(cell)) == "<class 'molsystem.cell.Cell'>"
    assert cell == [3.03, 3.03, 3.03, 90.0, 90.0, 90.0]

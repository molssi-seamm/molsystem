#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cell` in the `molsystem` package."""

import pytest  # noqa: F401


def test_construction(configuration):
    """Simplest test that we can make a Cell object"""
    configuration.periodicity = 3
    cell = configuration.cell
    assert str(type(cell)) == "<class 'molsystem.cell._Cell'>"


def test_construction_error(configuration):
    """Test that we get an error for nonperiodic systems"""
    with pytest.raises(TypeError) as e:
        configuration.cell
    if str(e.value) != "The configuration is not periodic!":
        print(e.value)
    assert str(e.value) == "The configuration is not periodic!"


def test_periodic_system(configuration):
    """Test making a simple periodic system."""
    configuration.periodicity = 3
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (3.03, 3.03, 3.03, 90, 90, 90)
    configuration.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol="V")
    assert configuration.atoms.n_atoms == 2
    assert configuration.version == 0


def test_periodic_context(configuration):
    """Test using context."""
    with configuration as tmp:
        tmp.periodicity = 3
        tmp.coordinate_system = "fractional"
        tmp.cell.parameters = (3.03, 3.03, 3.03, 90, 90, 90)
        tmp.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol="V")
    assert configuration.atoms.n_atoms == 2
    assert configuration.version == 1


def test_cell_object(vanadium):
    """Test getting the cell as a Cell object."""
    cell = vanadium.cell
    assert str(type(cell)) == "<class 'molsystem.cell._Cell'>"
    assert cell == [3.03, 3.03, 3.03, 90.0, 90.0, 90.0]

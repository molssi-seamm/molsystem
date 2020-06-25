#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cell` in the `molsystem` package."""

import pytest  # noqa: F401

from molsystem import Cell

pcell_parameters = [3, 4, 5, 60, 70, 80]


@pytest.fixture()
def cell():
    """A Cell object
    """
    cell = Cell()
    return cell


@pytest.fixture()
def pcell():
    """A periodic Cell object
    """
    cell = Cell()
    cell.periodicity = 3
    cell.parameters = pcell_parameters
    return cell


def test_construction():
    """Simplest test that we can make an Cell object"""
    system = None
    cell = Cell(system)
    assert str(type(cell)) == "<class 'molsystem.cell.Cell'>"


def test_str(cell):
    """String representation of the default cell."""
    assert str(cell) == 'not periodic'


def test_repr(cell):
    """representation of the default cell."""
    assert repr(cell) == (
        "{'periodicity': 0, 'cell': [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]}"
    )


def test_periodic_str(cell):
    """String representation of the default periodic cell."""
    cell.periodicity = 3
    assert str(cell) == (
        '3-D periodic, cell = 10.000 10.000 10.000 90.0 90.0 90.0'
    )


def test_version(cell):
    """The initial version of the cell"""
    assert cell.version == 0


def test_version_update(cell):
    """Version after an update."""
    with cell as tmp:
        tmp.periodicity = 3
    assert cell.version == 1


def test_get_cell_parameters(pcell):
    """Get the cell parameters as list."""
    assert pcell.parameters == pcell_parameters


def test_set_cell_parameters(cell):
    """Set the cell parameters using a list."""
    with cell as tmp:
        tmp.periodicity = 3
        tmp.parameters = pcell_parameters

    assert cell.parameters == pcell_parameters


def test_set_cell_parameters_scalar(cell):
    """Set the cell parameters using a single integer."""
    with pytest.raises(
        TypeError,
        match=('The cell parameters must be an iterable of floats, length 6.')
    ):
        with cell as tmp:
            tmp.periodicity = 3
            tmp.parameters = 1
    assert cell.version == 0


def test_set_cell_parameters_wrong_size(cell):
    """Set the cell parameters using a single integer."""
    with pytest.raises(
        TypeError,
        match=('The cell parameters must be an iterable of length 6.')
    ):
        with cell as tmp:
            tmp.periodicity = 3
            tmp.parameters = (9, 9, 9, 90)
    assert cell.version == 0


def test_set_cell_parameters_strings(cell):
    """Set the cell parameters using a list of strings."""
    with pytest.raises(
        TypeError, match='The cell parameters must be 6 floats.'
    ):
        with cell as tmp:
            tmp.periodicity = 3
            tmp.parameters = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')
    assert cell.version == 0


def test_set_cell_parameters_integers(cell):
    """Set the cell parameters using a list of integers."""
    with cell as tmp:
        tmp.periodicity = 3
        tmp.parameters = (10, 11, 12, 90, 90, 90)
    assert cell.version == 1


def test_get_a(pcell):
    """Get the a cell parameter."""
    assert pcell.a == pcell_parameters[0]


def test_set_a(pcell):
    """Set the a cell parameter."""
    with pcell as tmp:
        tmp.a = 20

    assert pcell.a == 20


def test_get_b(pcell):
    """Get the b cell parameter."""
    assert pcell.b == pcell_parameters[1]


def test_set_b(pcell):
    """Set the b cell parameter."""
    with pcell as tmp:
        tmp.b = 20

    assert pcell.b == 20


def test_get_c(pcell):
    """Get the c cell parameter."""
    assert pcell.c == pcell_parameters[2]


def test_set_c(pcell):
    """Set the c cell parameter."""
    with pcell as tmp:
        tmp.c = 20

    assert pcell.c == 20


def test_get_alpha(pcell):
    """Get the alpha cell parameter."""
    assert pcell.alpha == pcell_parameters[3]


def test_set_alpha(pcell):
    """Set the alpha cell parameter."""
    with pcell as tmp:
        tmp.alpha = 20

    assert pcell.alpha == 20


def test_get_beta(pcell):
    """Get the beta cell parameter."""
    assert pcell.beta == pcell_parameters[4]


def test_set_beta(pcell):
    """Set the beta cell parameter."""
    with pcell as tmp:
        tmp.beta = 20

    assert pcell.beta == 20


def test_get_gamma(pcell):
    """Get the gamma cell parameter."""
    assert pcell.gamma == pcell_parameters[5]


def test_set_gamma(pcell):
    """Set the gamma cell parameter."""
    with pcell as tmp:
        tmp.gamma = 20

    assert pcell.gamma == 20


def test_copy_constructor(pcell):
    """Get the gamma cell parameter."""
    cell2 = Cell(pcell)
    assert cell2 == pcell


def test_periodicity_error(cell):
    """Invalid periodicity value."""
    with pytest.raises(
        ValueError, match=r"The periodicity must be 0, 1, 2 or 3, not '4'"
    ):
        cell.periodicity = 4


def test_periodicity_error_type(cell):
    """Invalid periodicity value type."""
    with pytest.raises(
        TypeError, match=r"The periodicity must be an integer."
    ):
        cell.periodicity = 'abc'


def test_periodicity_error_notimplemented(cell):
    """Some periodicities not implemented yet."""
    with pytest.raises(
        NotImplementedError,
        match=r'1-D and 2-D periodicity not implemented yet.'
    ):
        cell.periodicity = 2


def test_no_change(pcell):
    """Update that does not change the cell."""
    with pcell as tmp:
        tmp.periodicity = 3

    assert pcell.version == 0

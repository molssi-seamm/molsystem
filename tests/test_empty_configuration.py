#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""An empty configuration's coordinates convert, and its periodicity changes (#121)."""

import numpy
import pytest  # noqa: F401


def test_empty_periodic_configuration_made_non_periodic(configuration):
    configuration.periodicity = 3
    configuration.cell.parameters = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    configuration.periodicity = 0
    assert configuration.periodicity == 0
    assert configuration.atoms.n_atoms == 0
    assert configuration.atoms.get_coordinates() == []


def test_empty_coordinates_convert(configuration):
    configuration.periodicity = 3
    cell = configuration.cell
    cell.parameters = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    assert cell.to_fractionals([]) == []
    assert cell.to_cartesians([]) == []
    assert cell.to_fractionals([], as_array=True).shape == (0, 3)
    assert cell.to_cartesians([], as_array=True).shape == (0, 3)
    # and non-empty ones as before
    uvw = cell.to_fractionals([[1.0, 2.0, 3.0]], as_array=True)
    assert uvw == pytest.approx(numpy.array([[0.1, 0.2, 0.3]]))

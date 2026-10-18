#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The element masses are seamm_util's, which are mendeleev's."""

import pytest
from seamm_util import element_data

from molsystem import elements


def test_one_table():
    assert elements.data is element_data
    for symbol, entry in element_data.items():
        assert elements.symbol_to_mass[symbol] == entry["atomic weight"]
        assert elements.atno_to_mass[entry["atomic number"]] == entry["atomic weight"]


@pytest.mark.parametrize(
    "symbol, mass", [("Li", 6.94), ("Ar", 39.948), ("F", 18.998403163)]
)
def test_masses(symbol, mass):
    assert elements.masses([symbol]) == [mass]


def test_argon_density():
    """1000 Ar in 47597.1 Å^3: 1.39368 g/mL (1.39122, 0.18% low, with the midpoint)."""
    mass = 1000 * elements.masses(["Ar"])[0]
    assert mass / 47597.1 * 1.0e24 / 6.02214076e23 == pytest.approx(1.39368, abs=1e-5)

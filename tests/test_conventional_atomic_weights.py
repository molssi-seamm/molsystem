#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Elements with an interval of atomic weights use CIAAW's conventional value."""

import pytest

CONVENTIONAL = {
    "H": 1.008,
    "Li": 6.94,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "Mg": 24.305,
    "Si": 28.085,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.95,
    "Br": 79.904,
    "Tl": 204.38,
    "Pb": 207.2,
}


@pytest.mark.parametrize("symbol", sorted(CONVENTIONAL))
def test_conventional_atomic_weight(symbol):
    from molsystem import elements

    assert elements.masses([symbol]) == [CONVENTIONAL[symbol]]


def test_argon_density():
    """1000 Ar in 47597.1 Å^3: 1.39375 g/mL (1.39122, 0.18% low, with the midpoint)."""
    from molsystem import elements

    mass = 1000 * elements.masses(["Ar"])[0]
    assert mass / 47597.1 * 1.0e24 / 6.02214076e23 == pytest.approx(1.39375, abs=1e-5)

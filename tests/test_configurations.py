#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for configurations."""

import pytest  # noqa: F401


def test_formula(AceticAcid):
    """Test the formula generation."""
    formula, empirical_formula, Z = AceticAcid.formula
    assert formula == "C2 H4 O2"
    assert empirical_formula == "C H2 O"
    assert Z == 2


def test_formula_periodic(vanadium):
    """Test the formula generation."""
    formula, empirical_formula, Z = vanadium.formula
    assert "".join(formula) == "V2"
    assert "".join(empirical_formula) == "V"
    assert Z == 2


def test_density(vanadium):
    """Test the density, and implicitly the mass and volume."""

    assert abs(vanadium.density - 6.0817308915133) < 1.0e-06

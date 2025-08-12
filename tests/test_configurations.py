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


def test_copy(H2O):
    """Test copying and changing a configuration."""
    original = H2O
    coordinates = [
        [0.0, 0.0, 0.0],
        [0.0, 0.8, -0.8],
        [0.0, -0.8, -0.8],
    ]
    copy = original.system.copy_configuration(original, name="Copy")

    assert original.coordinates == copy.coordinates

    copy.coordinates = coordinates

    assert original.coordinates != copy.coordinates


def test_copy_periodic(vanadium):
    """Test copying and changing a configuration."""
    original = vanadium

    copy = original.system.copy_configuration(original, name="Copy")

    assert original.cell.parameters == copy.cell.parameters

    copy.cell.parameters = [3.03, 3.03, 3.03, 90, 90, 91]

    assert original.cell.parameters != copy.cell.parameters


def test_spin(AceticAcid):
    """Test the spin multiplicity."""
    assert AceticAcid.spin_multiplicity == 1
    assert AceticAcid.spin_state == "singlet"
    AceticAcid.spin_state = "Triplet"
    assert AceticAcid.spin_multiplicity == 3

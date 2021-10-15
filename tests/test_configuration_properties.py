#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for properties."""

import pprint
import pytest  # noqa: F401


def test_empty(AceticAcid):
    """Test empty properties db."""
    assert AceticAcid.properties.known_properties() == []


def test_add(AceticAcid):
    """Test adding properties db."""
    properties = AceticAcid.properties
    assert not properties.exists("dipole moment")
    properties.add("dipole moment")
    assert properties.known_properties() == ["dipole moment"]
    assert properties.exists("dipole moment")


def test_put(AceticAcid):
    """Test putting some properties and values in."""
    properties = AceticAcid.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)
    result = properties.get("dipole moment")
    assert result == 2.0


def test_list(AceticAcid):
    """Test a simple query"""
    properties = AceticAcid.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)

    assert properties.list() == ["dipole moment"]


def test_get_all(AceticAcid):
    """Test getting all the properties and values."""
    answer = {
        "canonical SMILES": "CC(=O)O",
        "dipole moment": 2.0,
        "number of atoms": 8,
    }
    properties = AceticAcid.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)

    properties.add(
        "number of atoms", "int", description="The number of atoms in the system"
    )
    properties.put("number of atoms", AceticAcid.n_atoms)

    properties.add(
        "canonical SMILES", "str", description="The canonical SMILEs for the structure"
    )
    properties.put("canonical SMILES", AceticAcid.canonical_smiles)

    result = properties.get()
    if result != answer:
        pprint.pprint(result)
    assert result == answer

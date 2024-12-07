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
    assert result["dipole moment"]["value"] == 2.0


def test_list(AceticAcid):
    """Test a simple query"""
    properties = AceticAcid.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)

    assert properties.list() == ["dipole moment"]


def test_get_all(AceticAcid):
    """Test getting all the properties and values."""
    answer = {
        "canonical SMILES": {"cid": 1, "sid": 1, "value": "CC(=O)O"},
        "dipole moment": {"cid": 1, "sid": 1, "value": 2.0},
        "number of atoms": {"cid": 1, "sid": 1, "value": 8},
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


def test_list_property_ids(AceticAcid):
    """Test listing properties"""
    system = AceticAcid.system
    sys_properties = system.properties
    all_system_properties = []
    all_system_properties.append(sys_properties.add("dipole", "json", units="debye"))
    sys_properties.put("dipole", [0.0, 1.0, 2.0])
    all_system_properties.append(
        sys_properties.add("enthalpy of formation", "float", units="kJ/mol")
    )
    sys_properties.put("enthalpy of formation", -100.0)

    properties = AceticAcid.properties
    all_configuration_properties = []
    for code in ("Gaussian", "Psi4"):
        for basis in ("STO-3G", "6-31G"):
            all_configuration_properties.append(
                properties.add(f"dipole#{code}#HF/{basis}", "json", units="debye")
            )
            properties.put(f"dipole#{code}#HF/{basis}", [0.1, 0.9, 2.1])
            all_configuration_properties.append(
                properties.add(
                    f"enthalpy of formation#{code}#HF/{basis}", "float", units="kJ/mol"
                )
            )
            properties.put(f"enthalpy of formation#{code}#HF/{basis}", -100.1)

    assert sorted(sys_properties.list(as_ids=True)) == sorted(all_system_properties)
    assert sorted(properties.list(as_ids=True)) == sorted(all_configuration_properties)
    assert (
        sorted(properties.list(as_ids=True, include_system_properties=True))
    ) == sorted(all_configuration_properties + all_system_properties)

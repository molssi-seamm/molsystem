#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for properties."""

import pprint

import pytest  # noqa: F401


def test_empty(empty_db):
    """Test empty properties db."""
    assert empty_db.properties.known_properties() == []


def test_add(empty_db):
    """Test adding properties db."""
    properties = empty_db.properties
    assert not properties.exists("dipole moment")
    properties.add("dipole moment")
    assert properties.known_properties() == ["dipole moment"]
    assert properties.exists("dipole moment")


def test_put(AceticAcid):
    """Test putting some properties and values in."""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole moment")
    properties.put(cid, "dipole moment", 2.0)
    result = properties.get(cid, "dipole moment")
    assert result["dipole moment"]["value"] == 2.0


def test_query(AceticAcid):
    """Test a simple query"""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole moment")
    properties.put(cid, "dipole moment", 2.0)

    assert properties.query(what=["dipole moment"]) == {"dipole moment": [2.0]}


def test_duplicate_put(AceticAcid):
    """Test putting twice ... should only get one back with last value"""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole moment")
    properties.put(cid, "dipole moment", 2.0)
    properties.put(cid, "dipole moment", 3.0)

    assert properties.query(what=["dipole moment"]) == {"dipole moment": [3.0]}


def test_system_put(AceticAcid):
    """Test putting some properties and values in a system."""
    system = AceticAcid.system
    properties = system.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)
    result = properties.get("dipole moment")
    assert result["dipole moment"]["value"] == 2.0


def test_duplicate_system_put(AceticAcid):
    """Test putting some properties and values in a system twice."""
    system = AceticAcid.system
    properties = system.properties
    properties.add("dipole moment")
    properties.put("dipole moment", 2.0)
    properties.put("dipole moment", 3.0)
    result = properties.get("dipole moment")
    assert result["dipole moment"]["value"] == 3.0


def test_query2(properties):
    """Test a simple query"""
    result = properties.query(
        "float_0",
        "between",
        -10.0,
        -9.9,
        what=["float_8", "str_9"],
    )
    assert result["str_9"] == [
        "string 400",
        "string 401",
        "string 402",
        "string 403",
        "string 404",
        "string 405",
        "string 406",
        "string 407",
        "string 408",
        "string 409",
        "string 410",
    ]
    assert result["float_8"] == [
        -2.0,
        -1.99,
        -1.98,
        -1.97,
        -1.96,
        -1.95,
        -1.94,
        -1.93,
        -1.92,
        -1.91,
        -1.9,
    ]


def test_put_json(AceticAcid):
    """Test putting some properties and values in."""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole", "json")
    properties.put(cid, "dipole", [0.0, 1.0, 2.0])
    result = properties.get(cid, "dipole")
    assert result["dipole"]["value"] == [0.0, 1.0, 2.0]


def test_query_json(AceticAcid):
    """Test a simple query"""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole", "json")
    properties.put(cid, "dipole", [0.0, 1.0, 2.0])

    assert properties.query(what=["dipole"]) == {"dipole": [[0.0, 1.0, 2.0]]}


def test_list_properties(AceticAcid):
    """Test listing properties"""
    system = AceticAcid.system
    sys_properties = system.properties
    sys_properties.add("dipole", "json", units="debye")
    sys_properties.put("dipole", [0.0, 1.0, 2.0])
    sys_properties.add("enthalpy of formation", "float", units="kJ/mol")
    sys_properties.put("enthalpy of formation", -100.0)
    all_system_properties = ["dipole", "enthalpy of formation"]

    properties = AceticAcid.properties
    all_configuration_properties = []
    for code in ("Gaussian", "Psi4"):
        for basis in ("STO-3G", "6-31G"):
            properties.add(f"dipole#{code}#HF/{basis}", "json", units="debye")
            properties.put(f"dipole#{code}#HF/{basis}", [0.1, 0.9, 2.1])
            all_configuration_properties.append(f"dipole#{code}#HF/{basis}")
            properties.add(
                f"enthalpy of formation#{code}#HF/{basis}", "float", units="kJ/mol"
            )
            properties.put(f"enthalpy of formation#{code}#HF/{basis}", -100.1)
            all_configuration_properties.append(
                f"enthalpy of formation#{code}#HF/{basis}"
            )

    assert sorted(sys_properties.list()) == sorted(all_system_properties)
    assert sorted(properties.list()) == sorted(all_configuration_properties)
    assert sorted(properties.list(include_system_properties=True)) == sorted(
        all_configuration_properties + all_system_properties
    )

    correct = [
        "dipole#Gaussian#HF/STO-3G",
        "enthalpy of formation#Gaussian#HF/STO-3G",
        "dipole#Gaussian#HF/6-31G",
        "enthalpy of formation#Gaussian#HF/6-31G",
    ]

    if properties.list("*Gaussian*") != correct:
        pprint.pprint(properties.list("*Gaussian*"))
    assert properties.list("*Gaussian*") == correct


def test_get_properties(AceticAcid):
    """Test getting properties"""
    system = AceticAcid.system
    sys_properties = system.properties
    sys_properties.add("dipole", "json", units="debye")
    sys_properties.put("dipole", [0.0, 1.0, 2.0])
    sys_properties.add("enthalpy of formation", "float", units="kJ/mol")
    sys_properties.put("enthalpy of formation", -100.0)

    properties = AceticAcid.properties
    incr = 0.0
    for code in ("Gaussian", "Psi4"):
        for basis in ("STO-3G", "6-31G"):
            incr += 0.01
            properties.add(f"dipole#{code}#HF/{basis}", "json", units="debye")
            properties.put(
                f"dipole#{code}#HF/{basis}",
                [round(0.1 + incr, 3), round(0.9 - incr, 3), 2.1],
            )
            properties.add(
                f"enthalpy of formation#{code}#HF/{basis}", "float", units="kJ/mol"
            )
            properties.put(
                f"enthalpy of formation#{code}#HF/{basis}", round(-100.1 + incr, 3)
            )

    correct = {
        "dipole": {"cid": None, "sid": 1, "value": [0.0, 1.0, 2.0]},
        "enthalpy of formation": {"cid": None, "sid": 1, "value": -100.0},
    }

    if sys_properties.get() != correct:
        pprint.pprint(sys_properties.get())
    assert sys_properties.get() == correct

    correct = {
        "dipole#Gaussian#HF/6-31G": {"cid": 1, "sid": 1, "value": [0.12, 0.88, 2.1]},
        "dipole#Gaussian#HF/STO-3G": {"cid": 1, "sid": 1, "value": [0.11, 0.89, 2.1]},
        "dipole#Psi4#HF/6-31G": {"cid": 1, "sid": 1, "value": [0.14, 0.86, 2.1]},
        "dipole#Psi4#HF/STO-3G": {"cid": 1, "sid": 1, "value": [0.13, 0.87, 2.1]},
        "enthalpy of formation#Gaussian#HF/6-31G": {
            "cid": 1,
            "sid": 1,
            "value": -100.08,
        },
        "enthalpy of formation#Gaussian#HF/STO-3G": {
            "cid": 1,
            "sid": 1,
            "value": -100.09,
        },
        "enthalpy of formation#Psi4#HF/6-31G": {"cid": 1, "sid": 1, "value": -100.06},
        "enthalpy of formation#Psi4#HF/STO-3G": {"cid": 1, "sid": 1, "value": -100.07},
    }

    if properties.get() != correct:
        pprint.pprint(properties.get())
    assert properties.get() == correct

    correct = {
        "dipole": {"cid": None, "sid": 1, "value": [0.0, 1.0, 2.0]},
        "dipole#Gaussian#HF/6-31G": {"cid": 1, "sid": 1, "value": [0.12, 0.88, 2.1]},
        "dipole#Gaussian#HF/STO-3G": {"cid": 1, "sid": 1, "value": [0.11, 0.89, 2.1]},
        "dipole#Psi4#HF/6-31G": {"cid": 1, "sid": 1, "value": [0.14, 0.86, 2.1]},
        "dipole#Psi4#HF/STO-3G": {"cid": 1, "sid": 1, "value": [0.13, 0.87, 2.1]},
        "enthalpy of formation": {"cid": None, "sid": 1, "value": -100.0},
        "enthalpy of formation#Gaussian#HF/6-31G": {
            "cid": 1,
            "sid": 1,
            "value": -100.08,
        },
        "enthalpy of formation#Gaussian#HF/STO-3G": {
            "cid": 1,
            "sid": 1,
            "value": -100.09,
        },
        "enthalpy of formation#Psi4#HF/6-31G": {"cid": 1, "sid": 1, "value": -100.06},
        "enthalpy of formation#Psi4#HF/STO-3G": {"cid": 1, "sid": 1, "value": -100.07},
    }

    if properties.get(include_system_properties=True) != correct:
        pprint.pprint(properties.get(include_system_properties=True))
    assert properties.get(include_system_properties=True) == correct

    correct = {
        "dipole#Gaussian#HF/6-31G": {"cid": 1, "sid": 1, "value": [0.12, 0.88, 2.1]},
        "dipole#Gaussian#HF/STO-3G": {"cid": 1, "sid": 1, "value": [0.11, 0.89, 2.1]},
        "enthalpy of formation#Gaussian#HF/6-31G": {
            "cid": 1,
            "sid": 1,
            "value": -100.08,
        },
        "enthalpy of formation#Gaussian#HF/STO-3G": {
            "cid": 1,
            "sid": 1,
            "value": -100.09,
        },
    }

    if properties.get("*Gaussian*") != correct:
        pprint.pprint(properties.get("*Gaussian*"))
    assert properties.get("*Gaussian*") == correct

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for properties."""

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
    assert result == 2.0


def test_query(AceticAcid):
    """Test a simple query"""
    cid = AceticAcid.id
    properties = AceticAcid.system_db.properties
    properties.add("dipole moment")
    properties.put(cid, "dipole moment", 2.0)

    assert properties.query(what=["dipole moment"]) == {"dipole moment": [2.0]}


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

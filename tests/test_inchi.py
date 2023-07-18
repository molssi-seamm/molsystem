#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

"""Tests for handling InChI."""


def test_to_inchi(AceticAcid):
    """Create a InChI string from a system"""
    correct = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
    inchi = AceticAcid.to_inchi()

    if inchi != correct:
        print(inchi)
    assert inchi == correct


def test_to_inchikey(AceticAcid):
    """Create a InChIKey string from a system"""
    correct = "QTBSBXVTEAMEQO-UHFFFAOYSA-N"
    inchi = AceticAcid.to_inchi(key=True)

    if inchi != correct:
        print(inchi)
    assert inchi == correct


def test_to_inchi_with_name(AceticAcid):
    """Create a InChI string from a system, adding the name"""
    correct = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
    inchi = AceticAcid.inchi

    if inchi != correct:
        print(inchi)
    assert inchi == correct


def test_several_molecules(CH3COOH_3H2O):
    """System with acetic acid and 3 waters"""
    system = CH3COOH_3H2O
    correct = "InChI=1S/C2H4O2.3H2O/c1-2(3)4;;;/h1H3,(H,3,4);3*1H2"

    inchi = system.to_inchi()

    if inchi != correct:
        print(inchi)
    assert inchi == correct


def test_from_inchi(configuration):
    """Create a configuration from InChI"""
    correct = "CC(=O)O"
    configuration.from_inchi(
        "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)", name="acetic acid"
    )
    result = configuration.to_smiles(canonical=True)

    if result != correct:
        print(result)

    assert result == correct


def test_from_inchikey(configuration):
    """Create a configuration from an InChIKey"""
    correct = "CC(=O)O"

    configuration.from_inchikey("QTBSBXVTEAMEQO-UHFFFAOYSA-N")
    result = configuration.to_smiles(canonical=True)

    if result != correct:
        print(result)

    assert result == correct

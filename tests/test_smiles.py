#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

"""Tests for handling SMILES."""


def test_to_smiles(AceticAcid):
    """Create a SMILES string from a system"""
    correct = "CC(=O)O"
    smiles = AceticAcid.to_smiles()

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_to_smiles_with_name(AceticAcid):
    """Create a SMILES string from a system"""
    correct = "CC(=O)O"
    smiles = AceticAcid.smiles

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_to_canonical_smiles(AceticAcid):
    """Create a SMILES string from a system"""
    correct = "CC(=O)O"
    smiles = AceticAcid.to_smiles(canonical=True)

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_to_canonical_smiles_openbabel(AceticAcid):
    """Create a SMILES string from a system"""
    correct = "CC(=O)O"
    smiles = AceticAcid.to_smiles(canonical=True, flavor="openbabel")

    if smiles != correct:
        print(smiles)
    assert smiles == correct


@pytest.mark.openeye
def test_to_canonical_smiles_openeye(AceticAcid):
    """Create a SMILES string from a system"""
    correct = "CC(=O)O"
    smiles = AceticAcid.to_smiles(canonical=True, flavor="openeye")

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_several_molecules(CH3COOH_3H2O):
    """System with acetic acid and 3 waters"""
    system = CH3COOH_3H2O
    correct = "CC(=O)O.O.O.O"

    smiles = system.to_smiles()

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_several_molecules_canonical(CH3COOH_3H2O):
    """System with acetic acid and 3 waters"""
    system = CH3COOH_3H2O
    correct = "CC(=O)O.O.O.O"

    smiles = system.to_smiles(canonical=True)

    if smiles != correct:
        print(smiles)
    assert smiles == correct


def test_from_smiles(configuration):
    """Create a configuration from SMILES"""
    correct = "CC(=O)O"
    configuration.from_smiles("OC(=O)C", name="acetic acid")
    result = configuration.to_smiles(canonical=True)

    if result != correct:
        print(result)

    assert result == correct
    assert configuration.name == "acetic acid"


def test_from_smiles_openbabel(configuration):
    """Create a configuration from SMILES"""
    correct = "CC(=O)O"
    configuration.from_smiles("OC(=O)C", name="acetic acid", flavor="openbabel")
    result = configuration.to_smiles(canonical=True)

    if result != correct:
        print(result)

    assert result == correct
    assert configuration.name == "acetic acid"


@pytest.mark.openeye
def test_from_smiles_openeye(configuration):
    """Create a configuration from SMILES"""
    correct = "CC(=O)O"
    configuration.from_smiles("OC(=O)C", name="acetic acid", flavor="openeye")
    result = configuration.to_smiles(canonical=True)

    if result != correct:
        print(result)

    assert result == correct
    assert configuration.name == "acetic acid"

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


def test_rdkit_structure_independent_of_history(tmp_path):
    """The same SMILES gives the same structure whatever was built before it,
    in this process or another one."""
    import subprocess
    import sys
    import textwrap

    script = textwrap.dedent("""
        from molsystem import SystemDB

        db = SystemDB(filename="file:history?mode=memory&cache=shared")
        xyz = []
        for smiles in ("CCCCCCCC", "CCO", "c1ccccc1", "CCCCCCCC"):
            c = db.create_system(name=smiles).create_configuration(name=smiles)
            c.from_smiles(smiles, flavor="rdkit")
            if smiles == "CCCCCCCC":
                xyz.append(c.atoms.get_coordinates())
        assert xyz[0] == xyz[1], "differs within a process"
        print(xyz[0][0])
        """)
    first = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    ).stdout
    other_order = script.replace('("CCCCCCCC", "CCO",', '("CCO", "CCCCCCCC",')
    second = subprocess.run(
        [sys.executable, "-c", other_order], capture_output=True, text=True, check=True
    ).stdout
    assert first == second

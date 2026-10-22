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


@pytest.mark.parametrize("flavor", ["rdkit", "openbabel"])
@pytest.mark.parametrize("smiles", ["[Na+].O", "O.O.O", "[Li+].F[B-](F)(F)F"])
def test_fragments_do_not_overlap(configuration, flavor, smiles):
    """The fragments of a SMILES are placed apart: both toolkits embedded each at the
    origin, so they overlapped (2026-10-21: Na 0.38 Å from O in '[Na+].O')."""
    import numpy as np
    from molsystem.smiles import FRAGMENT_GAP

    configuration.from_smiles(smiles, flavor=flavor)
    xyz = np.array(configuration.atoms.coordinates)
    molecules = configuration.find_molecules(as_indices=True)
    assert len(molecules) == smiles.count(".") + 1
    for i, a in enumerate(molecules):
        for b in molecules[i + 1 :]:
            d = np.linalg.norm(xyz[a][:, None] - xyz[b][None], axis=2)
            assert d.min() >= FRAGMENT_GAP - 1.0e-6


def test_one_fragment_is_unchanged(configuration):
    """A single molecule is left where it was built."""
    import numpy as np

    configuration.from_smiles("CCO", reorient=False)
    before = np.array(configuration.atoms.coordinates)
    configuration._separate_fragments()
    np.testing.assert_allclose(configuration.atoms.coordinates, before)

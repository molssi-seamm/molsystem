#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `molsystem` package.

These tests center on angiotensin, a small peptide.
"""

from pathlib import Path
import pprint  # noqa: F401
import pytest  # noqa: F401

from molsystem import SystemDB

path = Path(__file__).resolve().parent
data_path = path / "data"

sequence = "H-Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-OH"

# SMILES for angiotensin II
# H-Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-OH
SMILES = (
    "[NH2][C@@H](CC(=O)O)C(=O)"
    "[NH][C@@H](CCCNC(=[NH2])N)C(=O)"
    "[NH][C@@H](C(C)C)C(=O)"
    "[NH][C@@H](Cc1ccc(cc1)O)C(=O)"
    "[NH][C@@H]([C@@H](CC)C)C(=O)"
    "[NH][C@@H](CC1=CNC=[NH]1)C(=O)"
    "[N]1[C@@H](CCC1)C(=O)"
    "[NH][C@@H](Cc1ccccc1)C(=O)O"
)
SMILES = (
    "[NH2][C@@H](CC(=O)O)C(=O)"
    "[NH][C@@H](CCCNC(=N)N)C(=O)"
    "[NH][C@@H](C(C)C)C(=O)"
    "[NH][C@@H](Cc1ccc(cc1)O)C(=O)"
    "[NH][C@@H]([C@@H](CC)C)C(=O)"
    "[NH][C@@H](CC1=CNC=N1)C(=O)"
    "[N]1[C@@H](CCC1)C(=O)"
    "[NH][C@@H](Cc1ccccc1)C(=O)O"
)
sidechains = {
    "ALA_LL": "[CH3X4]",
    "ARG_LL": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH2X3]",
    "ARG_LL_DHH12": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH1X2])[NH2X3]",
    "ARG_LL_DHH22": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH1X2]",
    "ARG_LL_RNH1": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH2X3])[NH2X3]",
    "ASN_LL": "[CH2X4][CH0X3](=[OH0X1])[NH2X3]",
    "ASP_LL": "[CH2X4][CH0X3](=[OH0X1])[OH1X2]",
    "ASP_LL_DHD2": "[CH2X4][CH0X3](=[OH0X1])[OH0X1]",
    "CYS_LL": "[CH2X4][SH1X2]",
    "CYS_LL_DHG": "[CH2X4][SH0X1]",
    "GLN_LL": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[NH2X3]",
    "GLU_LL": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH1X2]",
    "GLU_LL_DHE2": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH0X1]",
    "HIS_LL": "[CH2X4][cH0X3][nH1X3][cH1X3][nH1X3][cH1X3]",
    "HIS_LL_DHD1": "[CH2X4][cH0X3][nH0X2][cH1X3][nH1X3][cH1X3]",
    "HIS_LL_DHE2": "[CH2X4][cH0X3][nH1X3][cH1X3][nH0X2][cH1X3]",
    "ILE_LL": "[CH1X4]([CH2X4][CH3X4])[CH3X4]",
    "LEU_LL": "[CH2X4][CH1X4]([CH3X4])[CH3X4]",
    "LYS_LL": "[CH2X4][CH2X4][CH2X4][CH2X4][NH3X4]",
    "LYS_LL_DHZ3": "[CH2X4][CH2X4][CH2X4][CH2X4][NH2X3]",
    "MET_LL": "[CH2X4][CH2X4][SH0X2][CH3X4]",
    "PHE_LL": "[CH2X4][cH0X3][cH1X3][cH1X3][cH1X3][cH1X3][cH1X3]",
    "SER_LL": "[CH2X4][OH1X2]",
    "SER_LL_DHG": "[CH2X4][OH0X1]",
    "THR_LL": "[CH1X4]([OH1X2])[CH3X4]",
    "THR_LL_DHG1": "[CH1X4]([OH0X1])[CH3X4]",
    "TRP_LL": (
        "[CH2X4][CH0X3]=[CH1X3][NH1X3][CH0X3]=[CH0X3][CH1X3]=[CH1X3]" "[CH1X3]=[CH1X3]"
    ),
    "TRP_LL_DHE1": (
        "[CH2X4][CH0X3]=[CH1X3][NH0X2][CH0X3]=[CH0X3][CH1X3]=[CH1X3]" "[CH1X3]=[CH1X3]"
    ),
    "TYR_LL": "[CH2X4][cH0X3]1[cH1X3][cH1X3][cH0X3]([OH1X2])[cH1X3][cH1X3]1",
    "TYR_LL_DHH": ("[CH2X4][CH0X3]=[CH1X3][CH1X3]=[CH0X3]([CH1X3]=[CH1X3])[OH0X1]"),
    "VAL_LL": "[CH1X4]([CH3X4])[CH3X4]",
}

full_smarts = {
    "GLY_LL": "[NH1X2][CH2X4][CH0X2]=[OH0X1]",
    "PRO_LL": "[NH0X2][CH1X4]([CH0X2]=[OH0X1])[CH2X4][CH2X4][CH2X4]",
}


@pytest.fixture(scope="module")
def testdb():
    """Create a database with amino acids and templates for each."""
    db = SystemDB(filename="file:aa_templates_db?mode=memory&cache=shared")
    db.read_cif_file(data_path / "aa-variants-v1.cif")

    # Make the templates
    templates = db.templates
    for system in db.systems:
        tmp = system.name.split("_")
        if len(tmp) == 1:
            category = "amino acid"
        else:
            if tmp[1] == "LL":
                category = "residue"
            elif tmp[1] == "LEO2":
                category = "C-terminal residue"
            elif tmp[1] == "LEO2H":
                category = "protonated C-terminal residue"
            elif tmp[1] == "LFOH":
                category = "amino acid free neutral"
            elif tmp[1] == "LFZW":
                category = "amino acid free zwitterion"
            elif tmp[1] == "LSN3":
                category = "protonated N-terminal residue"
            else:
                raise ValueError(f"Don't recognize {tmp[1]} in {system.name}")

        templates.create(
            name=system.name, category=category, configuration=system.configuration.id
        )

    # Angiotensin II from SMILES
    system = db.create_system("angiotensin")
    configuration = system.create_configuration("SMILES")
    configuration.from_smiles(SMILES)

    # And read in
    system.read_cif_file(data_path / "1n9v.cif")

    yield db

    db.close()


@pytest.fixture()
def angiotensin(configuration):
    """A system with one configuration, containing angiotensin."""
    configuration.from_smiles(SMILES)

    return configuration


def test_correctness(angiotensin):
    """Check that the structure is correct."""
    assert angiotensin.formula[0] == "C50 H71 N13 O12"
    assert abs(angiotensin.mass - 1046.2) < 0.1


def test_testdb(testdb):
    """Test that the test database is correct."""
    db = testdb
    # There are 218 structures in aa-variants, plus the 1 angiotensin system
    assert db.n_systems == 219

    templates = db.templates
    assert templates.n_templates == 218

    system = db.get_system("angiotensin")
    # There is the SMILES configuration and 21 from the ensemble
    assert system.n_configurations == 22

    smiles = system.get_configuration("SMILES")

    assert smiles.atoms.n_atoms == 146
    assert smiles.formula[0] == "C50 H71 N13 O12"
    assert abs(smiles.mass - 1046.2) < 0.1

    angiotensin = system.get_configuration("representative")
    assert angiotensin.atoms.n_atoms == 146
    assert angiotensin.formula[0] == "C50 H71 N13 O12"
    assert abs(angiotensin.mass - 1046.2) < 0.1


def test_all_residue_search(testdb):
    """Testing locating residues in a peptide."""
    residues = {
        "ARG_LL_DHH12": [(9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)],
        "ARG_LL_DHH22": [(9, 10, 11, 12, 13, 14, 15, 17, 16, 18, 19)],
        "HIS_LL_DHD1": [(47, 48, 49, 50, 54, 53, 52, 51, 55, 56)],
        "ILE_LL": [(39, 40, 41, 42, 43, 44, 45, 46)],
        "PHE_LL": [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74)],
        "TYR_LL": [(27, 28, 29, 30, 31, 32, 33, 36, 34, 35, 37, 38)],
        "VAL_LL": [(20, 21, 22, 23, 24, 25, 26)],
    }
    n_terminal = {
        "ASP_LL": [(1, 2, 3, 4, 5, 6, 7, 8)],
    }
    c_terminal = {
        "PHE_LL": [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75)],
    }

    db = testdb

    if False:
        print("SMARTS")
        for name in db.templates.names("residue"):
            tmp = db.templates.get(name, "residue")
            print(f"'{name}': '{tmp.smarts}',")

        print("SMILES")
        for name in db.templates.names("residue"):
            tmp = db.templates.get(name, "residue")
            print(f"'{name}': '{tmp.smiles}',")

    system = db.get_system("angiotensin")
    angiotensin = system.get_configuration("SMILES")

    # print(angiotensin.smiles)

    for name, sc in sidechains.items():
        smarts = f"[NH1X3][C@@H]({sc})[CX3]=[OX1]"
        result = angiotensin.find_substructures(smarts)
        if len(result) == 0:
            assert name not in residues
        else:
            if name not in residues or result != residues[name]:
                print(f"'{name}': {result},")
            assert result == residues[name]

    for name, sc in sidechains.items():
        smarts = f"[NH2][C@@H]({sc})[CX3]=[OX1]"
        result = angiotensin.find_substructures(smarts)
        if len(result) == 0:
            assert name not in n_terminal
        else:
            if name not in n_terminal or result != n_terminal[name]:
                print(f"'{name}': {result},")
            assert result == n_terminal[name]

    for name, sc in sidechains.items():
        smarts = f"[NX3][C@@H]({sc})[CX3](=[OX1])[OH1X2]"
        result = angiotensin.find_substructures(smarts)
        if len(result) == 0:
            assert name not in c_terminal
        else:
            if name not in c_terminal or result != c_terminal[name]:
                print(f"'{name}': {result},")
            assert result == c_terminal[name]

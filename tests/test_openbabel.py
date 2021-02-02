#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenBabel mixin of the class."""

import pprint  # noqa: F401

import pytest  # noqa: F401

# H-Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-OH
SMILES = (
    '[NH2][C@@H](CC(=O)O)C(=O)'
    '[NH][C@@H](CCCNC(=[NH2])N)C(=O)'
    '[NH][C@@H](C(C)C)C(=O)'
    '[NH][C@@H](Cc1ccc(cc1)O)C(=O)'
    '[NH][C@@H]([C@@H](CC)C)C(=O)'
    '[NH][C@@H](CC1=CNC=[NH]1)C(=O)'
    '[N]1[C@@H](CCC1)C(=O)'
    '[NH][C@@H](Cc1ccccc1)C(=O)O'
)
sidechains = {
    'ALA_LL': '[CH3X4]',
    'ARG_LL': '[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH2X3]',
    'ARG_LL_DHH12': '[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH1X2])[NH2X3]',
    'ARG_LL_DHH22': '[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH1X2]',
    'ARG_LL_RNH1': '[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH2X3])[NH2X3]',
    'ASN_LL': '[CH2X4][CH0X3](=[OH0X1])[NH2X3]',
    'ASP_LL': '[CH2X4][CH0X3](=[OH0X1])[OH1X2]',
    'ASP_LL_DHD2': '[CH2X4][CH0X3](=[OH0X1])[OH0X1]',
    'CYS_LL': '[CH2X4][SH1X2]',
    'CYS_LL_DHG': '[CH2X4][SH0X1]',
    'GLN_LL': '[CH2X4][CH2X4][CH0X3](=[OH0X1])[NH2X3]',
    'GLU_LL': '[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH1X2]',
    'GLU_LL_DHE2': '[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH0X1]',
    'HIS_LL': '[CH2X4][CH0X3][NH1X3]=[CH1X3][NH1X3][CH1X3]',
    'HIS_LL_DHD1': '[CH2X4][CH0X3][NH0X2]=[CH1X3][NH1X3][CH1X3]',
    'HIS_LL_DHE2': '[CH2X4][CH0X3][NH1X3]=[CH1X3][NH0X2][CH1X3]',
    'ILE_LL': '[CH1X4]([CH2X4][CH3X4])[CH3X4]',
    'LEU_LL': '[CH2X4][CH1X4]([CH3X4])[CH3X4]',
    'LYS_LL': '[CH2X4][CH2X4][CH2X4][CH2X4][NH3X4]',
    'LYS_LL_DHZ3': '[CH2X4][CH2X4][CH2X4][CH2X4][NH2X3]',
    'MET_LL': '[CH2X4][CH2X4][SH0X2][CH3X4]',
    'PHE_LL': '[CH2X4][cH0X3][cH1X3][cH1X3][cH1X3][cH1X3][cH1X3]',
    'SER_LL': '[CH2X4][OH1X2]',
    'SER_LL_DHG': '[CH2X4][OH0X1]',
    'THR_LL': '[CH1X4]([OH1X2])[CH3X4]',
    'THR_LL_DHG1': '[CH1X4]([OH0X1])[CH3X4]',
    'TRP_LL':
        (
            '[CH2X4][CH0X3]=[CH1X3][NH1X3][CH0X3]=[CH0X3][CH1X3]=[CH1X3]'
            '[CH1X3]=[CH1X3]'
        ),
    'TRP_LL_DHE1':
        (
            '[CH2X4][CH0X3]=[CH1X3][NH0X2][CH0X3]=[CH0X3][CH1X3]=[CH1X3]'
            '[CH1X3]=[CH1X3]'
        ),
    'TYR_LL': '[CH2X4][cH0X3]1[cH1X3][cH1X3][cH0X3]([OH1X2])[cH1X3][cH1X3]1',
    'TYR_LL_DHH':
        ('[CH2X4][CH0X3]=[CH1X3][CH1X3]=[CH0X3]([CH1X3]=[CH1X3])[OH0X1]'),
    'VAL_LL': '[CH1X4]([CH3X4])[CH3X4]'
}

full_smarts = {
    'GLY_LL': '[NH1X2][CH2X4][CH0X2]=[OH0X1]',
    'PRO_LL': '[NH0X2][CH1X4]([CH0X2]=[OH0X1])[CH2X4][CH2X4][CH2X4]'
}


def test_substructure(CH3COOH_3H2O):
    """Test the finding substructures in the configuration."""
    answer1 = [(5, 6, 7)]
    answer2 = [(5, 6, 7, 8)]
    answer3 = [(9,), (12,), (15,)]
    answer4 = [(10, 9, 11), (13, 12, 14), (16, 15, 17)]
    answer5 = [(9, 10, 11), (12, 13, 14), (15, 16, 17)]

    configuration = CH3COOH_3H2O

    # Just the C and O of the carboxyl, not to H on O
    result = configuration.find_substructures('C(=O)O')
    if result != answer1:
        pprint.pprint(result)
    assert result == answer1

    # All four atoms of the carboxyl group
    result = configuration.find_substructures('C(=O)[O][H]')
    if result != answer2:
        pprint.pprint(result)
    assert result == answer2

    # The Oxygen atoms of the waters
    result = configuration.find_substructures('[OH2]')
    if result != answer3:
        pprint.pprint(result)
    assert result == answer3

    # All the atoms in waters, order H-O-H
    result = configuration.find_substructures('[H][O][H]')
    if result != answer4:
        pprint.pprint(result)
    assert result == answer4

    # All the atoms in waters, order O-H-H
    result = configuration.find_substructures('[O]([H])[H]')
    if result != answer5:
        pprint.pprint(result)
    assert result == answer5


def test_substructure_ordering(disordered):
    """Test the ordering of atoms in subsets."""
    answer1 = 'CC(=O)O'
    answer2 = [(1, 5, 6, 7), (16, 12, 11, 10)]
    answer3 = 'C([H])([H])([H])C(=O)O[H]'
    answer4 = [(1, 2, 3, 4, 5, 6, 7, 8), (16, 15, 14, 13, 12, 11, 10, 9)]

    configuration = disordered
    templates = configuration.create_molecule_templates(create_subsets=False)

    # Without hydrogens
    smiles = templates[0].smiles
    if smiles != answer1:
        pprint.pprint(smiles)
    assert smiles == answer1

    result = configuration.find_substructures(smiles)
    if result != answer2:
        pprint.pprint(result)
    assert result == answer2

    # With hydrogens
    smiles = templates[0].to_smiles(hydrogens=True)
    if smiles != answer3:
        pprint.pprint(smiles)
    assert smiles == answer3

    result = configuration.find_substructures(smiles)
    if result != answer4:
        pprint.pprint(result)
    assert result == answer4


def test_all_residue_search(configuration):
    """Testing locating residues in a peptide."""
    residues = {
        'ARG_LL': [(9, 10, 11, 12, 13, 14, 15, 17, 16, 18, 19)],
        'ARG_LL_RNH1': [(9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)],
        'HIS_LL': [(47, 48, 49, 50, 54, 53, 52, 51, 55, 56)],
        'ILE_LL': [(39, 40, 41, 42, 43, 44, 45, 46)],
        'PHE_LL': [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74)],
        'TYR_LL': [(27, 28, 29, 30, 31, 32, 33, 36, 34, 35, 37, 38)],
        'VAL_LL': [(20, 21, 22, 23, 24, 25, 26)],
    }
    n_terminal = {
        'ASP_LL': [(1, 2, 3, 4, 5, 6, 7, 8)],
    }
    c_terminal = {
        'PHE_LL': [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75)],
    }

    configuration.from_smiles(SMILES)

    for name, sc in sidechains.items():
        smarts = f"[NH1X3][C@@H]({sc})[CX3]=[OX1]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in residues
        else:
            if result != residues[name]:
                print(f"'{name}': {result},")
            assert result == residues[name]

    for name, sc in sidechains.items():
        smarts = f"[NH2][C@@H]({sc})[CX3]=[OX1]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in n_terminal
        else:
            if result != n_terminal[name]:
                print(f"'{name}': {result},")
            assert result == n_terminal[name]

    for name, sc in sidechains.items():
        smarts = f"[NX3][C@@H]({sc})[CX3](=[OX1])[OH1X2]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in c_terminal
        else:
            if result != c_terminal[name]:
                print(f"'{name}': {result},")
            assert result == c_terminal[name]

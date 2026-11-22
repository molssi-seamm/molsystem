#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Symmetry-aware RMSD must not swap the two O atoms of an acid group: two
rotamers of a COOH are different conformers (Thermochemistry, 2026-11-21)."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import molsystem

ALA_ALA = "C[C@H](N)C(=O)N[C@@H](C)C(=O)O"


def rotate(points, a, b, angle):
    """Rotate points by angle about the axis from a to b."""
    k = (b - a) / np.linalg.norm(b - a)
    p = points - a
    c, s = np.cos(angle), np.sin(angle)
    return a + p * c + np.cross(k, p) * s + np.outer(p @ k, k) * (1 - c)


def acid_rotamers():
    """Ala-Ala, and the same with its COOH turned 180 degrees about the Ca-C bond."""
    mol = Chem.AddHs(Chem.MolFromSmiles(ALA_ALA))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    xyz = mol.GetConformer().GetPositions()
    ca, c, o1, o2, h = mol.GetSubstructMatch(
        Chem.MolFromSmarts("[CX4][CX3](=O)[OX2][H]")
    )
    flipped = xyz.copy()
    moved = [o1, o2, h]
    flipped[moved] = rotate(xyz[moved], xyz[ca], xyz[c], np.pi)
    other = Chem.RWMol(mol)
    conformer = other.GetConformer()
    for i, point in enumerate(flipped):
        conformer.SetAtomPosition(i, point.tolist())
    return Chem.RWMol(mol), other


def test_acid_rotamers_differ_without_h():
    a, b = acid_rotamers()
    result = molsystem.RMSD(b, a, symmetry=True)
    assert result["RMSD"] > 0.5  # 0.95 A; RDKit's default gave 0.10


def test_acid_rotamers_differ_with_h():
    a, b = acid_rotamers()
    result = molsystem.RMSD(b, a, symmetry=True, include_h=True)
    assert result["RMSD"] > 0.5


def test_true_symmetry_still_found():
    """Swapping the numbering of a methyl group's H atoms is still recognized."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    xyz = mol.GetConformer().GetPositions()
    hs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(0).GetNeighbors() if n.GetSymbol() == "H"
    ]
    swapped = xyz.copy()
    swapped[[hs[0], hs[1]]] = xyz[[hs[1], hs[0]]]
    other = Chem.RWMol(mol)
    conformer = other.GetConformer()
    for i, point in enumerate(swapped):
        conformer.SetAtomPosition(i, point.tolist())
    result = molsystem.RMSD(other, Chem.RWMol(mol), symmetry=True, include_h=True)
    assert result["RMSD"] < 1.0e-6


def test_rdkit_option_can_be_turned_on():
    """RDKit's terminal-group symmetrization is still available on request."""
    a, b = acid_rotamers()
    result = molsystem.RMSD(b, a, symmetry=True, symmetrize_terminal_groups=True)
    assert result["RMSD"] < 0.5

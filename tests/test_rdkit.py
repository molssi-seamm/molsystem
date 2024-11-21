#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the RDKit mixin of the class."""

import pprint  # noqa: F401
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
import pytest  # noqa: F401


def test_to_RDKMol(Acetate):
    """Test creating a RDKMol object from a structure."""
    correct = {
        "float property": 3.14,
        "float property,units": "kcal/mol",
        "int property": 2,
        "net charge": -1,
        "spin multiplicity": 1,
        "str property": "Hi!",
    }

    mol = Acetate.to_RDKMol(properties="all")

    bond_types = {
        rdkit.Chem.rdchem.BondType.SINGLE: 1,
        rdkit.Chem.rdchem.BondType.DOUBLE: 2,
        rdkit.Chem.rdchem.BondType.TRIPLE: 3,
    }
    bond_list = [bond_types[bt.GetBondType()] for bt in mol.GetBonds()]

    assert Acetate.n_atoms == mol.GetNumAtoms()
    assert Acetate.bonds.bondorders == bond_list

    data = mol.GetPropsAsDict()
    if data != correct:
        pprint.pprint(data)
    assert data == correct


def test_from_RDKMol(configuration):
    """Test creating a structure from an RDKMol object."""
    mol = rdkit.Chem.MolFromSmiles("CC=O")
    mol = rdkit.Chem.rdmolops.AddHs(mol, addCoords=True)

    rdkit.Chem.AllChem.EmbedMolecule(mol)

    configuration.from_RDKMol(mol)

    assert configuration.n_atoms == 7
    assert configuration.bonds.bondorders == [1, 2, 1, 1, 1, 1]

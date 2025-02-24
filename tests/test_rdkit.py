#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the RDKit mixin of the class."""

import pprint  # noqa: F401
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
import pytest  # noqa: F401

from molsystem import rdkit_version


def test_version():
    """Test the version of RDKit."""
    rdkit_version()


def test_to_RDKMol(Acetate):
    """Test creating a RDKMol object from a structure."""
    correct = {
        "SEAMM|XYZ|json|": "[\n"
        "    [1.0797, 0.0181, -0.0184],\n"
        "    [0.5782, 3.1376, 0.2813],\n"
        "    [0.7209, -0.6736, -0.7859],\n"
        "    [0.7052, -0.3143, 0.9529],\n"
        "    [0.5713, 1.3899, -0.3161],\n"
        "    [-0.1323, 1.7142, -1.2568],\n"
        "    [0.9757, 2.297, 0.5919]\n"
        "]",
        "SEAMM|configuration name|str|": "acetate",
        "SEAMM|float property|float|kcal/mol": 3.14,
        "SEAMM|int property|int|": 2,
        "SEAMM|net charge|int|": -1,
        "SEAMM|spin multiplicity|int|": 1,
        "SEAMM|str property|str|": "Hi!",
        "SEAMM|system name|str|": "acetate",
    }

    mol = Acetate.to_RDKMol(properties="*")

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

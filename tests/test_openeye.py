#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenEye mixin of the class."""

import pprint  # noqa: F401
import pytest  # noqa: F401

try:
    from openeye import oechem, oeomega
else:
    pass


@pytest.mark.openeye
def test_to_OEGraphMol(configuration):
    """Test creating a OEGraphMol object from a structure."""
    mol = configuration.to_OEGraphMol()

    assert configuration.n_atoms == mol.NumAtoms()
    bond_list = [bond.GetOrder() for bond in mol.GetBonds()]
    assert configuration.bonds.bondorders == bond_list


@pytest.mark.openeye
def test_from_OEMol(configuration):
    """Test creating a structure from an RDKMol object."""
    mol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol, "CC=O")

    omega = oeomega.OEOmega()
    # Only generate one conformer for our molecule
    omega.SetMaxConfs(1)
    # Set to False to pick random stereoisomer if stereochemistry is not specified
    omega.SetStrictStereo(False)
    # Be a little loose about atom typing to ensure parameters are available to omega
    # for all molecules
    omega.SetStrictAtomTypes(False)

    # Add hydrogens
    oechem.OEAddExplicitHydrogens(mol, False, True)

    assert mol.NumAtoms() == 7

    configuration.from_OEMol(mol)

    assert configuration.n_atoms == 7
    assert configuration.bonds.bondorders == [1, 2, 1, 1, 1, 1]

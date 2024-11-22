#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenEye mixin of the class."""

import pprint  # noqa: F401
import pytest  # noqa: F401

try:
    from openeye import oechem, oeomega
except ImportError:
    pass

from molsystem import openeye_version


@pytest.mark.openeye
def test_version():
    """Test the version of the OpenEye toolkit."""
    openeye_version()


@pytest.mark.openeye
def test_to_OEGraphMol(Acetate):
    """Test creating a OEGraphMol object from a structure."""
    correct = {
        "float property": 3.14,
        "float property,units": "kcal/mol",
        "int property": 2,
        "net charge": -1,
        "spin multiplicity": 1,
        "str property": "Hi!",
    }

    mol = Acetate.to_OEGraphMol(properties="all")

    assert Acetate.n_atoms == mol.NumAtoms()
    bond_list = [bond.GetOrder() for bond in mol.GetBonds()]
    assert Acetate.bonds.bondorders == bond_list

    data = {}
    for tmp in mol.GetDataIter():
        tag = tmp.GetTag()
        attribute = oechem.OEGetTag(tag)
        value = mol.GetData(tag)
        data[attribute] = value
    if data != correct:
        pprint.pprint(data)
    assert data == correct


@pytest.mark.openeye
def test_from_OEMol(configuration):
    """Test creating a structure from an RDKMol object."""
    correct = {"float property": 3.14, "integer property": 42, "string property": "foo"}

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

    # Add some properties
    mol.SetIntData("integer property", 42)
    mol.SetDoubleData("float property", 3.14)
    mol.SetStringData("string property", "foo")
    mol.SetIntData("net charge", -2)
    mol.SetIntData("spin multiplicity", 3)

    configuration.from_OEMol(mol)

    assert configuration.n_atoms == 7
    assert configuration.bonds.bondorders == [1, 2, 1, 1, 1, 1]

    result = configuration.properties.get()
    if result != correct:
        pprint.pprint(result)
    assert result == correct

    assert configuration.charge == -2
    assert configuration.spin_multiplicity == 3

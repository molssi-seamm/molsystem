#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the SystemDB class.

These tests are standalone, and test that the class heirarchy works to
a minimal extent, and that we can make a system with atoms and bonds.

Once these test work then conftest.py should function and the other
tests can explore the functionlaity in more detail.
"""

import math
import pytest  # noqa: F401

from molsystem import SystemDB


def test_construction():
    """Simplest test that we can make a SystemDB object"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    assert str(type(db)) == "<class 'molsystem.system_db.SystemDB'>"
    assert db.n_systems == 0
    db.close()


def test_system():
    """Test that we can get a System object"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    assert str(type(system)) == "<class 'molsystem.system._System'>"
    assert system.n_configurations == 0
    db.close()


def test_configuration():
    """Test that we can get a Configuration object"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration()
    _type = str(type(configuration))
    assert _type == "<class 'molsystem.configuration._Configuration'>"
    db.close()


def test_atoms():
    """Test that we can get an Atoms object"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration()
    atoms = configuration.atoms
    assert str(type(atoms)) == "<class 'molsystem.atoms._Atoms'>"
    db.close()


def test_adding_atoms():
    """Test that we can add atoms"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration()
    atoms = configuration.atoms
    symbol = ["Ar"] * 4
    X = [2.0, 0.0, -2.0, 0.0]
    Y = [0.0, 2.0, 0.0, -2.0]
    Z = [0.0] * 4
    ids = atoms.append(symbol=symbol, x=X, y=Y, z=Z)
    assert atoms.n_atoms == 4
    assert ids == [1, 2, 3, 4]
    db.close()


def test_bonds():
    """Test that we can get an Bonds object"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration()
    bonds = configuration.bonds
    assert str(type(bonds)) == "<class 'molsystem.bonds._Bonds'>"
    db.close()


def test_adding_bonds():
    """Test that we can add bonds"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system(name="default")
    configuration = system.create_configuration()

    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    name = ["O", "H1", "H2"]
    i_atom = [0, 0]
    j_atom = [1, 2]

    atoms = configuration.atoms
    atoms.add_attribute("name", "str", default="")
    atom_ids = atoms.append(atno=atno, x=X, y=Y, z=Z, name=name)

    i = [atom_ids[x] for x in i_atom]
    j = [atom_ids[x] for x in j_atom]

    bonds = configuration.bonds
    bond_ids = bonds.append(i=j, j=i)  # flipped on purpose so code orders.

    assert bonds.n_bonds == 2
    assert bond_ids == [1, 2]

    db.close()


def test_adding_reading_cif_file(amino_acids):
    """Test that we can read a cif file with many systems."""
    answer = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLU",
        "GLN",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
    assert amino_acids.n_systems == 20
    names = amino_acids.names
    if names != answer:
        print(names)
    assert names == answer


def test_accessing_systems(amino_acids):
    """Test that we can access each system."""
    answer = [
        13,
        27,
        17,
        16,
        14,
        19,
        20,
        10,
        21,
        22,
        22,
        25,
        20,
        23,
        17,
        14,
        17,
        27,
        24,
        19,
    ]
    assert amino_acids.n_systems == 20
    n_atoms = [x.configuration.n_atoms for x in amino_acids.systems]
    if n_atoms != answer:
        print(n_atoms)
    assert n_atoms == answer

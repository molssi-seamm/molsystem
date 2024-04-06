#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `symmetry` in the `molsystem` package."""

import json
from pathlib import Path

import pytest  # noqa: F401
from molsystem.symmetry import _Symmetry

path = Path(__file__).resolve().parent
data_path = path / "data"


def test_construction(configuration):
    """Simplest test that we can make a Symmetry object"""
    symmetry = _Symmetry(configuration)
    assert str(type(symmetry)) == "<class 'molsystem.symmetry._Symmetry'>"


def test_spacegroup_names(symmetry):
    """Test the correspondance of spacegroup names and Hall numbers."""
    path = Path(__file__).parent / "data" / "spacegroup_names.json"
    with path.open() as fd:
        correct = json.load(fd)
    answer = symmetry.spacegroup_names_to_hall()
    if answer != correct:
        print(json.dumps(answer, indent=3))
    assert answer == correct


def test_spacegroup(symmetry):
    """Handling of a spacegroup."""
    symmetry.configuration.periodicity = 3
    symmetry.group = "P b c a"

    correct = [
        "x,y,z",
        "-x,-y,-z",
        "-x+1/2,-y,z+1/2",
        "x+1/2,y,-z+1/2",
        "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,z",
        "-x,y+1/2,-z+1/2",
        "x,-y+1/2,z+1/2",
    ]

    symops = symmetry.symops
    if symops != correct:
        print(symops)

    assert symops == correct
    assert symmetry.n_symops == len(correct)


def test_operations(symmetry):
    """The symmetry operations as matrices."""
    symmetry.configuration.periodicity = 3
    symmetry.group = "P b c a"

    correct = """\
[[[ 1.   0.   0.   0. ]
  [ 0.   1.   0.   0. ]
  [ 0.   0.   1.   0. ]
  [ 0.   0.   0.   1. ]]

 [[-1.   0.   0.   0. ]
  [ 0.  -1.   0.   0. ]
  [ 0.   0.  -1.   0. ]
  [ 0.   0.   0.   1. ]]

 [[-1.   0.   0.   0.5]
  [ 0.  -1.   0.   0. ]
  [ 0.   0.   1.   0.5]
  [ 0.   0.   0.   1. ]]

 [[ 1.   0.   0.   0.5]
  [ 0.   1.   0.   0. ]
  [ 0.   0.  -1.   0.5]
  [ 0.   0.   0.   1. ]]

 [[ 1.   0.   0.   0.5]
  [ 0.  -1.   0.   0.5]
  [ 0.   0.  -1.   0. ]
  [ 0.   0.   0.   1. ]]

 [[-1.   0.   0.   0.5]
  [ 0.   1.   0.   0.5]
  [ 0.   0.   1.   0. ]
  [ 0.   0.   0.   1. ]]

 [[-1.   0.   0.   0. ]
  [ 0.   1.   0.   0.5]
  [ 0.   0.  -1.   0.5]
  [ 0.   0.   0.   1. ]]

 [[ 1.   0.   0.   0. ]
  [ 0.  -1.   0.   0.5]
  [ 0.   0.   1.   0.5]
  [ 0.   0.   0.   1. ]]]"""

    operations = symmetry.symmetry_matrices
    if str(operations) != correct:
        print(operations)

    assert str(operations) == correct


def test_coordinates(diamond):
    """Test the creation of the symmetric atom coordinates."""
    correct = """\
[
    [
        0.0,
        0.0,
        0.0
    ],
    [
        0.0,
        0.5,
        0.5
    ],
    [
        0.25,
        0.25,
        0.25
    ],
    [
        0.25,
        0.75,
        0.75
    ],
    [
        0.5,
        0.0,
        0.5
    ],
    [
        0.5,
        0.5,
        0.0
    ],
    [
        0.75,
        0.25,
        0.75
    ],
    [
        0.75,
        0.75,
        0.25
    ]
]"""
    assert diamond.n_atoms == 8
    uvw = json.dumps(diamond.coordinates, indent=4)
    if uvw != correct:
        print(uvw)
    assert uvw == correct


def test_bonds(diamond):
    """Test the creation of the symmetric bonds."""
    correct = [1.5459] * 16
    # diamond.symmetry.loglevel = "DEBUG"
    # diamond.bonds.loglevel = "DEBUG"
    R = diamond.bonds.get_lengths(as_array=True)
    assert diamond.n_asymmetric_bonds == 4
    assert diamond.n_bonds == 16
    assert all(R.round(4) == correct)


def test_trigonal(configuration):
    """Test symmetry for a trigonal systems with x-y."""
    path = data_path / "trigonal.cif"
    cif_text = path.read_text()
    configuration.from_cif_text(cif_text)

    assert configuration.n_asymmetric_atoms == 95
    assert configuration.n_atoms == 3 * 95
    assert configuration.n_asymmetric_bonds == 100
    assert configuration.n_bonds == 3 * 100


def test_molecules(benzene):
    """Test getting the molecules benzene."""
    correct_indices = [
        [0, 7, 8, 15, 16, 23, 24, 31, 33, 38, 41, 46],
        [1, 6, 9, 14, 17, 22, 25, 30, 32, 39, 40, 47],
        [2, 5, 10, 13, 18, 21, 26, 29, 35, 36, 43, 44],
        [3, 4, 11, 12, 19, 20, 27, 28, 34, 37, 42, 45],
    ]

    # benzene.symmetry.loglevel = "DEBUG"
    # benzene.bonds.loglevel = "DEBUG"

    indices = benzene.find_molecules(as_indices=True)

    if indices != correct_indices:
        print(indices)

    assert indices == correct_indices


def test_inverse_operations(benzene):
    """Test getting the inverse of the symmetry operations."""
    correct = [0, 1, 2, 3, 4, 5, 6, 7]
    inverse_ops = benzene.symmetry.inverse_operations
    if inverse_ops != correct:
        print(inverse_ops)
    assert inverse_ops == correct


def test_symmetrizing_coordinates(benzene):
    """Test symmetrizing the coordinates to get asymmetric atom coordinates."""
    correct = [
        [0.1297, 0.5762, 0.40803],
        [0.2172, 0.6275, 0.346],
        [0.1235, 0.6328, 0.54518],
        [0.2068, 0.7225, 0.5756],
        [0.0057, 0.4432, 0.36289],
        [0.0095, 0.4051, 0.2704],
    ]
    atoms = benzene.atoms
    uvw = atoms.coordinates
    uvw_asym, delta = benzene.symmetry.symmetrize_coordinates(uvw)
    if uvw_asym != correct:
        print(uvw_asym)
    assert uvw_asym == correct

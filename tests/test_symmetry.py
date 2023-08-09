#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `symmetry` in the `molsystem` package."""

import json
from pathlib import Path

import pytest  # noqa: F401
from molsystem.symmetry import _Symmetry


def test_construction(configuration):
    """Simplest test that we can make a Symmetry object"""
    symmetry = _Symmetry(configuration)
    assert str(type(symmetry)) == "<class 'molsystem.symmetry._Symmetry'>"


def test_spacegroup_names(symmetry):
    """Test the correspondance of spacegroup names and Hall numbers."""
    path = Path(__file__).parent / "data" / "spacegroup_names.json"
    with path.open() as fd:
        correct = json.load(fd)
    answer = symmetry.spacegroup_names_to_hall
    if answer != correct:
        print(json.dumps(answer, indent=3))
    assert answer == correct


def test_spacegroup(symmetry):
    """Handling of a spacegroup."""
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
    assert symmetry.n_operations == len(correct)


def test_operations(symmetry):
    """The symmetry operations as matrices."""
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

    operations = symmetry.operations
    if str(operations) != correct:
        print(operations)

    assert str(operations) == correct

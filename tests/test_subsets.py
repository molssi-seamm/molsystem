#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of subsets."""

import pprint  # noqa: F401

import pytest  # noqa: F401


def test_construction(CH3COOH_3H2O):
    """Simplest test that we can make a new subsets object"""
    configuration = CH3COOH_3H2O
    subsets = configuration.subsets
    assert str(type(subsets)) == "<class 'molsystem.subsets._Subsets'>"


def test_n_subsets(CH3COOH_3H2O):
    """Test getting the number of subsets."""
    configuration = CH3COOH_3H2O
    subsets = configuration.subsets
    assert subsets.n_subsets == 0

    db = configuration.system_db
    templates = db.templates
    tid = templates.create(name="H2O", category="molecule")

    subset = subsets.create(template=tid)
    assert subsets.n_subsets == 1
    assert subset.id == 1

    ids = subsets.get_ids(tid)
    assert ids == [1]


def test_simple_subsets(simple_templates):
    """Test adding atoms to a simple subset."""
    answer = [6, 1, 1, 1, 6, 8, 8, 1]

    db = simple_templates
    templates = db.templates
    acy = templates.get("acy", "molecule")
    hoh = templates.get("hoh", "molecule")

    system = db.get_system("acetic acid")
    configuration = system.configuration
    subsets = configuration.subsets

    subsets.create(acy, atoms=[*range(1, 9)])
    subset = subsets.get(acy)[0]
    assert subset.atoms.n_atoms == 8

    assert subset.atoms.atomic_numbers == answer

    subsets.create(hoh, atoms=[9, 10, 11])
    subsets.create(hoh, atoms=[12, 13, 14])
    # for this one add the atoms later
    last = subsets.create(hoh)

    waters = subsets.get(hoh)

    # adding atoms to last water
    last.atoms.add([15, 16, 17])

    for water in waters:
        assert water.atoms.atomic_numbers == [8, 1, 1]
        assert water.atoms.symbols == ["O", "H", "H"]

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the topology mixin of the system class."""

import pprint

import pytest  # noqa: F401


def test_bonded_neighbors(AceticAcid):
    """Test the generation of bonded neighbors."""
    result = {
        1: [2, 3, 4, 5],
        2: [1],
        3: [1],
        4: [1],
        5: [1, 6, 7],
        6: [5],
        7: [5, 8],
        8: [7],
    }

    neighbors = AceticAcid.bonded_neighbors()
    if neighbors != result:
        pprint.pprint(neighbors)
    assert neighbors == result


def test_molecules(AceticAcid):
    """Test the finding molecules in the system."""
    result = [[1, 2, 3, 4, 5, 6, 7, 8]]

    molecules = AceticAcid.find_molecules()
    if molecules != result:
        pprint.pprint(molecules)
    assert molecules == result


def test_multiple_molecules(CH3COOH_3H2O):
    """Test the finding molecules in the system."""
    result = [[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11], [12, 13, 14], [15, 16, 17]]

    molecules = CH3COOH_3H2O.find_molecules()
    if molecules != result:
        pprint.pprint(molecules)
    assert molecules == result


def test_molecule_templates(disordered):
    """Test making templates for the molecules."""
    result_tids = [1]

    configuration = disordered
    templates = configuration.create_molecule_templates(create_subsets=False)
    tids = [x.id for x in templates]
    if tids != result_tids:
        print("tids")
        pprint.pprint(tids)
    assert tids == result_tids


def test_molecule_subsets(disordered):
    """Test making templates and subsets for the molecules."""
    result_tids = [1]
    result_sids = {1: [1, 2]}

    configuration = disordered
    templates, subsets = configuration.create_molecule_templates()
    tids = [x.id for x in templates]
    if tids != result_tids:
        print("tids")
        pprint.pprint(tids)
    sids = {t: [x.id for x in s] for t, s in subsets.items()}
    if sids != result_sids:
        print("sids")
        pprint.pprint(sids)
    assert tids == result_tids
    assert sids == result_sids

    # And that the atoms are indeed ordered correctly
    tid = tids[0]
    subset1, subset2 = subsets[tid]
    atnos1 = subset1.atoms.atomic_numbers
    atnos2 = subset2.atoms.atomic_numbers
    if atnos1 != atnos2:
        print(f"atnos1 = {atnos1}")
        print(f"atnos2 = {atnos2}")
    assert atnos1 == atnos2

    coords1 = subset1.atoms.coordinates
    coords2 = subset2.atoms.coordinates

    # At the moment the methyl hydrogens are not properly ordered -- they are
    # equivalent -- so remove them to test.
    c1 = list(coords1)
    del c1[1:4]
    c2 = list(coords2)
    del c2[1:4]

    if c1 != c2:
        print("coords1")
        pprint.pprint(coords1)
        print("coords2")
        pprint.pprint(coords2)
    assert c1 == c2

    assert subset1.atoms.get_n_atoms("atno", "==", 6) == 2

    count = 0
    for row in subset2.atoms.atoms("atno", "==", 6):
        count += 1
    assert count == 2

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
        8: [7]
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
    result = [
        [1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11], [12, 13, 14], [15, 16, 17]
    ]

    molecules = CH3COOH_3H2O.find_molecules()
    if molecules != result:
        pprint.pprint(molecules)
    assert molecules == result


def test_molecule_subsets(CH3COOH_3H2O):
    """Test making subsets for the molecules."""
    result = [2, 3, 4, 5]

    system = CH3COOH_3H2O
    sids = system.create_molecule_subsets()
    if sids != result:
        pprint.pprint(sids)
    assert sids == result

    sid = sids[1]
    assert system.atoms.n_atoms(subset=sid) == 3
    assert system.bonds.n_bonds(subset=sid) == 0  # No bonds defined!


def test_molecule_templates(disordered):
    """Test making templates for the molecules."""
    result_tids = [2, 2]
    result_sids = {2: [2, 3]}

    system = disordered
    tids, sids = system.create_molecule_templates()
    if tids != result_tids:
        print('tids')
        pprint.pprint(tids)
    if sids != result_sids:
        print('sids')
        pprint.pprint(sids)
    assert tids == result_tids
    assert sids == result_sids

    # And that the atoms are indeed order correctly
    tid = tids[0]
    sid1, sid2 = sids[tid]
    atnos1 = system.atoms.atomic_numbers(subset=sid1)
    atnos2 = system.atoms.atomic_numbers(subset=sid2, template_order=True)
    if atnos1 != atnos2:
        print(f'atnos1 = {atnos1}')
        print(f'atnos2 = {atnos2}')
    assert atnos1 == atnos2

    coords1 = system.atoms.coordinates(subset=sid1)
    coords2 = system.atoms.coordinates(subset=sid2, template_order=True)

    # At the moment the methyl hydrogens are not properly ordered -- they are
    # equivalent -- so remove them to test.
    c1 = list(coords1)
    del c1[1:4]
    c2 = list(coords2)
    del c2[1:4]

    if c1 != c2:
        print('coords1')
        pprint.pprint(coords1)
        print('coords2')
        pprint.pprint(coords2)
    assert c1 == c2

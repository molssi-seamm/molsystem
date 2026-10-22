#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Hypervalent and coordination centres from SMILES get their VSEPR geometry.

RDKit and Open Babel embedded PF6- with fluorines 0.32-1.19 Å apart (2026-10-22: a
Packmol box of LiPF6 died at the first MD step).
"""

import numpy as np
import pytest

from molsystem import vsepr

FLAVORS = ["rdkit", "openbabel"]


def _geometry(configuration):
    """The centre (most bonds), its bond lengths and ligand angles, and the closest
    pair of fluorines."""
    xyz = np.array(configuration.atoms.coordinates)
    symbols = configuration.atoms.symbols
    index = {a: k for k, a in enumerate(configuration.atoms.ids)}
    neighbors = {k: [] for k in range(len(symbols))}
    for i, j in zip(
        configuration.bonds.get_column_data("i"),
        configuration.bonds.get_column_data("j"),
    ):
        neighbors[index[i]].append(index[j])
        neighbors[index[j]].append(index[i])
    centre = max(neighbors, key=lambda k: (symbols[k] != "H", len(neighbors[k])))
    vectors = xyz[neighbors[centre]] - xyz[centre]
    lengths = np.linalg.norm(vectors, axis=1)
    unit = vectors / lengths[:, None]
    cosines = unit @ unit.T
    angles = np.degrees(
        np.arccos(np.clip(cosines[np.triu_indices(len(unit), 1)], -1, 1))
    )
    fluorines = [k for k, s in enumerate(symbols) if s == "F"]
    ff = [
        np.linalg.norm(xyz[a] - xyz[b])
        for i, a in enumerate(fluorines)
        for b in fluorines[i + 1 :]
    ]
    return symbols[centre], lengths, np.sort(angles), min(ff) if ff else None


@pytest.mark.parametrize("flavor", FLAVORS)
@pytest.mark.parametrize(
    "smiles, centre, length",
    [
        ("F[P-](F)(F)(F)(F)F", "P", 1.606),
        ("FS(F)(F)(F)(F)F", "S", 1.56),
        ("F[Si-2](F)(F)(F)(F)F", "Si", 1.68),
        ("F[As-](F)(F)(F)(F)F", "As", 1.72),
    ],
)
def test_octahedral(configuration, flavor, smiles, centre, length):
    configuration.from_smiles(smiles, flavor=flavor)
    symbol, lengths, angles, ff = _geometry(configuration)
    assert symbol == centre
    np.testing.assert_allclose(lengths, length, atol=1e-3)
    # 12 cis pairs at 90 and 3 trans pairs at 180
    np.testing.assert_allclose(angles[:12], 90.0, atol=0.5)
    np.testing.assert_allclose(angles[12:], 180.0, atol=0.5)
    assert ff >= vsepr.MIN_F_F


@pytest.mark.parametrize("flavor", FLAVORS)
@pytest.mark.parametrize("smiles", ["F[B-](F)(F)F", "[O-]Cl(=O)(=O)=O"])
def test_tetrahedral_is_left_alone(configuration, flavor, smiles):
    """Four electron pairs (BF4-, ClO4-) are embedded well, and left as they are."""
    configuration.from_smiles(smiles, flavor=flavor)
    _, _, angles, _ = _geometry(configuration)
    np.testing.assert_allclose(angles, 109.47, atol=3.0)


@pytest.mark.parametrize("flavor", FLAVORS)
def test_lone_pairs(configuration, flavor):
    """SF4 is a seesaw (lone pair equatorial), XeF4 square planar."""
    configuration.from_smiles("FS(F)(F)F", flavor=flavor)
    _, _, angles, _ = _geometry(configuration)
    np.testing.assert_allclose(angles, [90, 90, 90, 90, 120, 180], atol=0.5)
    configuration.from_smiles("F[Xe](F)(F)F", flavor=flavor)
    _, _, angles, _ = _geometry(configuration)
    np.testing.assert_allclose(angles, [90, 90, 90, 90, 180, 180], atol=0.5)


@pytest.mark.parametrize("flavor", FLAVORS)
def test_metal_complex_and_large_ligands(configuration, flavor):
    """[Co(NH3)6]3+ is octahedral; the methoxy groups of P(OMe)6- are turned to
    keep clear of each other (the contact check would otherwise raise)."""
    configuration.from_smiles("N[Co+3](N)(N)(N)(N)N", flavor=flavor)
    symbol, _, angles, _ = _geometry(configuration)
    assert symbol == "Co"
    np.testing.assert_allclose(angles[:12], 90.0, atol=0.5)
    configuration.from_smiles("CO[P-](OC)(OC)(OC)(OC)OC", flavor=flavor)
    symbol, _, angles, _ = _geometry(configuration)
    assert symbol == "P"
    np.testing.assert_allclose(angles[12:], 180.0, atol=0.5)


@pytest.mark.parametrize("flavor", FLAVORS)
def test_salt_of_two_fragments(configuration, flavor):
    """LiPF6: an octahedral anion, and the cation clear of it."""
    configuration.from_smiles("[Li+].F[P-](F)(F)(F)(F)F", flavor=flavor)
    _, lengths, angles, ff = _geometry(configuration)
    np.testing.assert_allclose(lengths, 1.606, atol=1e-3)
    assert ff >= vsepr.MIN_F_F
    assert configuration.charge == 0


def test_directions():
    """The polyhedra, with the lone pairs where VSEPR puts them."""
    seesaw = vsepr.ligand_directions(4, 1)
    assert sum(abs(v[2]) > 0.99 for v in seesaw) == 2  # both axial sites kept
    planar = vsepr.ligand_directions(4, 2)
    assert np.allclose(planar[:, 2], 0.0)  # the two lone pairs are trans
    assert np.allclose(vsepr.ligand_directions(4, square_planar=True)[:, 2], 0.0)


def test_contact_check():
    """Too-close atoms raise, naming them."""
    symbols = ["F", "P", "F"]
    xyz = [[0, 0, 0], [1.6, 0, 0], [0.3, 0, 0]]
    with pytest.raises(RuntimeError, match="F1-F3 0.30"):
        vsepr.check_contacts(symbols, xyz, [(0, 1, 1), (1, 2, 1)])

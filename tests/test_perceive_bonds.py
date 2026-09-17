# -*- coding: utf-8 -*-

"""Tests of TopologyMixin.perceive_bonds."""

import math

import numpy as np
import pytest

r_oh = 0.9572
half = math.radians(104.52 / 2)
XH = r_oh * math.sin(half)
ZH = r_oh * math.cos(half)


def _water(origin, rot=np.eye(3)):
    o = np.asarray(origin, dtype=float)
    return [o, o + rot @ np.array([XH, 0, ZH]), o + rot @ np.array([-XH, 0, ZH])]


def _add(configuration, coords, atnos, cell=None):
    xyz = np.array(coords)
    if cell is not None:
        configuration.periodicity = 3
        configuration.cell.parameters = cell
        configuration.coordinate_system = "Cartesian"
    return configuration.atoms.append(
        atno=atnos, x=xyz[:, 0].tolist(), y=xyz[:, 1].tolist(), z=xyz[:, 2].tolist()
    )


def test_water_dimer_nonperiodic(configuration):
    coords = _water([0, 0, 0]) + _water([2.9, 0, 0])
    _add(configuration, coords, [8, 1, 1, 8, 1, 1])
    n = configuration.perceive_bonds()
    assert n == 4
    assert configuration.bonds.n_bonds == 4
    mols = configuration.find_molecules(as_indices=True)
    assert mols == [[0, 1, 2], [3, 4, 5]]
    assert set(configuration.bonds.bondorders) == {1}


def test_existing_bonds_and_replace(configuration):
    _add(configuration, _water([0, 0, 0]), [8, 1, 1])
    assert configuration.perceive_bonds() == 2
    with pytest.raises(RuntimeError):
        configuration.perceive_bonds()
    assert configuration.perceive_bonds(replace=True) == 2
    assert configuration.bonds.n_bonds == 2


def test_periodic_molecule_across_boundary(configuration):
    """A water sitting on the cell boundary is bonded through the boundary."""
    L = 10.0
    # oxygen just inside the far face, hydrogens wrapped to the near face
    coords = _water(
        [L - 0.1, 5.0, 5.0], rot=np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    )
    coords = [c % L for c in coords]
    coords += _water([5.0, 5.0, 5.0])
    _add(configuration, coords, [8, 1, 1, 8, 1, 1], cell=[L, L, L, 90, 90, 90])
    assert configuration.perceive_bonds() == 4
    assert configuration.find_molecules(as_indices=True) == [[0, 1, 2], [3, 4, 5]]


def test_general_cell_matches_orthorhombic(configuration, two_configurations):
    """The supercell (general cell) path gives the same bonds as the KD-tree path."""
    rng = np.random.default_rng(3)
    L = 12.0
    coords, atnos = [], []
    for _ in range(40):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        a, b, c, d = q
        R = np.array(
            [
                [
                    a * a + b * b - c * c - d * d,
                    2 * (b * c - a * d),
                    2 * (b * d + a * c),
                ],
                [
                    2 * (b * c + a * d),
                    a * a - b * b + c * c - d * d,
                    2 * (c * d - a * b),
                ],
                [
                    2 * (b * d - a * c),
                    2 * (c * d + a * b),
                    a * a - b * b - c * c + d * d,
                ],
            ]
        )
        coords += _water(rng.uniform(0, L, 3), rot=R)
        atnos += [8, 1, 1]
    _add(configuration, coords, atnos, cell=[L, L, L, 90, 90, 90])
    n1 = configuration.perceive_bonds()
    b1 = sorted(
        tuple(sorted(p))
        for p in zip(*[configuration.bonds.get_column_data(k) for k in "ij"])
    )
    # Same geometry, cell declared marginally non-orthorhombic -> general path
    other = two_configurations[0]
    _add(other, coords, atnos, cell=[L, L, L, 90, 90, 90.0001])
    n2 = other.perceive_bonds()
    b2 = sorted(
        tuple(sorted(p)) for p in zip(*[other.bonds.get_column_data(k) for k in "ij"])
    )
    assert n1 == n2 > 0
    ids1 = {aid: k for k, aid in enumerate(configuration.atoms.ids)}
    ids2 = {aid: k for k, aid in enumerate(other.atoms.ids)}
    assert [(ids1[i], ids1[j]) for i, j in b1] == [(ids2[i], ids2[j]) for i, j in b2]


def test_ions_not_bonded_and_bf4(configuration):
    """Li+ stays a lone atom; BF4- gets exactly four B-F bonds, no F-F."""
    b = np.array([0.0, 0.0, 0.0])
    d = 1.40 / math.sqrt(3)
    fs = [
        b + d * np.array(v) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])
    ]
    li = [np.array([2.0, 0.0, 0.0])]  # ~1.9 A from the nearest F: would bond by radii
    _add(configuration, [b] + fs + li, [5, 9, 9, 9, 9, 3])
    assert configuration.perceive_bonds() == 4
    mols = configuration.find_molecules(as_indices=True)
    assert mols == [[0, 1, 2, 3, 4], [5]]
    # Opting in to bonding ions attaches the Li to the fluorines within reach
    n = configuration.perceive_bonds(replace=True, exclude=())
    assert n > 4
    assert configuration.find_molecules(as_indices=True) == [[0, 1, 2, 3, 4, 5]]


def test_hydrogen_limit_keeps_shortest(configuration):
    """A hydrogen between two oxygens bonds only to the nearer one."""
    coords = [[0.0, 0, 0], [1.0, 0, 0], [2.1, 0, 0]]  # O, H, O : both within tolerance
    _add(configuration, coords, [8, 1, 8])
    assert configuration.perceive_bonds() == 1
    i, j = configuration.bonds.get_column_data(
        "i"
    ), configuration.bonds.get_column_data("j")
    ids = configuration.atoms.ids
    assert sorted((ids.index(i[0]), ids.index(j[0]))) == [0, 1]
    # Raising the limit lets it bond to both
    assert configuration.perceive_bonds(replace=True, max_bonds={"H": 2}) == 2


def test_tolerance_and_radii_override(configuration):
    coords = [[0.0, 0, 0], [1.3, 0, 0]]  # O-H stretched to 1.3 A
    _add(configuration, coords, [8, 1])
    # default: 1.2 * (0.63 + 0.32) = 1.14 A < 1.3 A, so no bond
    assert configuration.perceive_bonds() == 0
    # a looser tolerance, or a larger radius for H, picks it up
    assert configuration.perceive_bonds(tolerance=1.4) == 1
    assert configuration.perceive_bonds(replace=True, radii={"H": 0.5}) == 1


def test_errors(configuration):
    _add(configuration, _water([0, 0, 0]), [8, 1, 1])
    with pytest.raises(ValueError):
        configuration.perceive_bonds(method="voronoi")


def _diamond_conventional(configuration, a=3.567):
    f = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
    f = f + [(x + 0.25, y + 0.25, z + 0.25) for x, y, z in f]
    configuration.periodicity = 3
    configuration.cell.parameters = [a, a, a, 90, 90, 90]
    configuration.coordinate_system = "fractional"
    configuration.atoms.append(
        atno=[6] * 8, x=[p[0] for p in f], y=[p[1] for p in f], z=[p[2] for p in f]
    )


def test_diamond_conventional_cell(configuration):
    """8-atom diamond: 16 bonds of 1.545 Å, some through the cell boundary."""
    _diamond_conventional(configuration)
    assert configuration.perceive_bonds() == 16
    lengths = np.array(configuration.bonds.get_lengths())
    assert len(lengths) == 16
    assert np.allclose(lengths, 3.567 * math.sqrt(3) / 4, atol=1e-3)
    # every carbon has four bonds
    counts = np.zeros(8, dtype=int)
    ids = list(configuration.atoms.ids)
    for i, j in zip(*[configuration.bonds.get_column_data(k) for k in "ij"]):
        counts[ids.index(i)] += 1
        counts[ids.index(j)] += 1
    assert counts.tolist() == [4] * 8
    assert configuration.find_molecules(as_indices=True) == [list(range(8))]


def test_diamond_primitive_cell(configuration):
    """2-atom diamond: four bonds between the same two atoms via different images."""
    a = 3.567 / math.sqrt(2)
    configuration.periodicity = 3
    configuration.cell.parameters = [a, a, a, 60, 60, 60]
    configuration.coordinate_system = "fractional"
    configuration.atoms.append(atno=[6, 6], x=[0, 0.25], y=[0, 0.25], z=[0, 0.25])
    assert configuration.perceive_bonds() == 4
    lengths = np.array(configuration.bonds.get_lengths())
    assert np.allclose(lengths, 3.567 * math.sqrt(3) / 4, atol=1e-3)
    symops = configuration.bonds.get_column_data("symop2")
    assert len(set(symops)) == 4  # four distinct images


def test_one_atom_cell_self_image_bonds(configuration):
    """Simple-cubic polonium: each atom bonds to its own images, 3 unique bonds."""
    a = 3.35
    configuration.periodicity = 3
    configuration.cell.parameters = [a, a, a, 90, 90, 90]
    configuration.coordinate_system = "fractional"
    configuration.atoms.append(atno=[84], x=[0.0], y=[0.0], z=[0.0])
    assert configuration.perceive_bonds() == 3
    assert np.allclose(configuration.bonds.get_lengths(), a, atol=1e-6)


def test_offsets_with_unwrapped_coordinates(configuration):
    """Stored coordinates outside [0,1) still give correct offsets and lengths."""
    L = 10.0
    coords = _water([L + 0.3, 5.0, 5.0])  # whole molecule sits past the far face
    coords[1] = coords[1] - np.array([L, 0, 0])  # one H stored in the home cell
    _add(configuration, coords, [8, 1, 1], cell=[L, L, L, 90, 90, 90])
    assert configuration.perceive_bonds() == 2
    assert np.allclose(configuration.bonds.get_lengths(), r_oh, atol=1e-6)

# -*- coding: utf-8 -*-

"""Sensible geometries for hypervalent and coordination centres built from SMILES.

RDKit and Open Babel embed centres with five or more electron pairs badly: PF6- from
``F[P-](F)(F)(F)(F)F`` had two fluorines 0.32-1.19 Å apart (2026-10-22: it killed an
MD run at its first step). This places the ligands of such centres on the ideal VSEPR
polyhedron, with the lone pairs where VSEPR puts them, at standard bond lengths, and
moves each ligand's other atoms with it.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

#: Bond lengths, Å, where the sum of covalent radii is poor
BOND_LENGTHS = {
    ("F", "P"): 1.606,
    ("B", "F"): 1.40,
    ("Al", "F"): 1.81,
    ("F", "Si"): 1.68,
    ("F", "S"): 1.56,
    ("As", "F"): 1.72,
    ("F", "Sb"): 1.88,
    ("Cl", "O"): 1.45,
}

#: Transition-metal ions with four ligands that are square planar (d8)
SQUARE_PLANAR = {
    ("Ni", 2),
    ("Pd", 2),
    ("Pt", 2),
    ("Au", 3),
    ("Rh", 1),
    ("Ir", 1),
}

#: Transition metals, whose ligands are placed by their number alone
TRANSITION_METALS = set(
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "Hf Ta W Re Os Ir Pt Au Hg".split()
)
#: Main-group atoms that can have five or more electron pairs with fewer than five
#: ligands, thanks to lone pairs (SF4, ClF3, XeF4, I3-, ...)
HYPERVALENT = set("P As Sb S Se Te Cl Br I Xe Kr".split())

#: The least distances, Å, between atoms more than two bonds apart: heavy atoms,
#: and any pair with a hydrogen; and between fluorines however far apart
MIN_HEAVY = 1.8
MIN_WITH_H = 1.2
MIN_F_F = 2.2

_T = 1.0 / np.sqrt(3.0)
_POLYHEDRA = {
    2: [(0, 0, 1), (0, 0, -1)],
    3: [(1, 0, 0), (-0.5, np.sqrt(3) / 2, 0), (-0.5, -np.sqrt(3) / 2, 0)],
    4: [(_T, _T, _T), (_T, -_T, -_T), (-_T, _T, -_T), (-_T, -_T, _T)],
    # trigonal bipyramid: two axial, then three equatorial
    5: [
        (0, 0, 1),
        (0, 0, -1),
        (1, 0, 0),
        (-0.5, np.sqrt(3) / 2, 0),
        (-0.5, -np.sqrt(3) / 2, 0),
    ],
    # octahedron, in trans pairs
    6: [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)],
    # pentagonal bipyramid: two axial, then five equatorial
    7: [(0, 0, 1), (0, 0, -1)]
    + [(np.cos(2 * np.pi * k / 5), np.sin(2 * np.pi * k / 5), 0) for k in range(5)],
}
#: Where the lone pairs go, in order, for each steric number: equatorial in a
#: trigonal bipyramid, trans in an octahedron, axial in a pentagonal bipyramid
_LONE_PAIR_SITES = {5: [2, 3, 4], 6: [0, 1, 2, 3], 7: [0, 1]}


def ligand_directions(n_ligands, n_lone_pairs=0, square_planar=False):
    """Unit vectors for the ligands of a centre, from VSEPR.

    Parameters
    ----------
    n_ligands : int
        The number of atoms bonded to the centre.
    n_lone_pairs : int
        Its lone pairs.
    square_planar : bool
        For a d8 metal with four ligands: an octahedron less a trans pair.

    Returns
    -------
    numpy.ndarray(n_ligands, 3)
    """
    if square_planar and n_ligands == 4:
        return np.array(_POLYHEDRA[6][2:], dtype=float)
    steric = n_ligands + n_lone_pairs
    if steric not in _POLYHEDRA:
        raise ValueError(f"No polyhedron for a steric number of {steric}.")
    sites = list(range(steric))
    for site in _LONE_PAIR_SITES.get(steric, list(range(steric)))[:n_lone_pairs]:
        sites.remove(site)
    return np.array([_POLYHEDRA[steric][i] for i in sites], dtype=float)


def _valence_electrons(group):
    """The valence electrons of a main-group element, from its group; None for a
    transition metal or an element without a group."""
    if group is None:
        return None
    if group <= 2:
        return group
    if group >= 13:
        return group - 10
    return None


def _rotation(a, b):
    """The rotation matrix taking unit vector a onto unit vector b."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1.0e-8:
        if c > 0:
            return np.eye(3)
        # 180°: about any axis perpendicular to a
        axis = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(axis) < 1.0e-8:
            axis = np.cross(a, [0.0, 1.0, 0.0])
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1.0 + c)


def _branch(neighbors, start, centre):
    """The atoms reachable from ``start`` without passing through ``centre``."""
    seen = {start}
    todo = [start]
    while todo:
        atom = todo.pop()
        for other in neighbors[atom]:
            if other != centre and other not in seen:
                seen.add(other)
                todo.append(other)
    return seen


def regularize(symbols, xyz, bonds, charges=None):
    """Place the ligands of hypervalent and coordination centres by VSEPR.

    A centre is reshaped when it has five or more electron pairs (a main-group atom
    with a steric number of 5 to 7: PF6-, SF6, SF4, ClF3, XeF4, IF7, ...) or is a
    transition metal with two or more ligands. Centres with four pairs or fewer
    (BF4-, ClO4-, ...) are left as the embedder made them. Chelating ligands are not
    handled: such centres are left alone, with a warning.

    Parameters
    ----------
    symbols : [str]
        The element of each atom.
    xyz : numpy.ndarray(n, 3)
        The coordinates, Å.
    bonds : [(int, int, int)]
        The bonds as (i, j, bond order), 0-based.
    charges : [int]
        The formal charge of each atom, or None for none.

    Returns
    -------
    (numpy.ndarray(n, 3), [int])
        The coordinates, and the atoms whose ligands were placed (empty if none).
    """
    n = len(symbols)
    xyz = np.array(xyz, dtype=float)
    charges = [0] * n if charges is None else [int(q or 0) for q in charges]
    neighbors = [[] for _ in range(n)]
    order_sum = [0.0] * n
    for i, j, order in bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)
        order = 1.5 if order == 5 else float(order)
        order_sum[i] += order
        order_sum[j] += order

    # Cheap screen first, so most molecules need no element data
    candidates = [
        i
        for i in range(n)
        if len(neighbors[i]) >= 5
        or (len(neighbors[i]) >= 2 and symbols[i] in TRANSITION_METALS | HYPERVALENT)
    ]
    if not candidates:
        return xyz, []

    import mendeleev

    elements = {}

    def element(symbol):
        if symbol not in elements:
            elements[symbol] = mendeleev.element(symbol)
        return elements[symbol]

    def bond_length(a, b):
        key = tuple(sorted((a, b)))
        if key in BOND_LENGTHS:
            return BOND_LENGTHS[key]
        return (
            element(a).covalent_radius_pyykko + element(b).covalent_radius_pyykko
        ) / 100.0

    changed = []
    for centre in candidates:
        ligands = neighbors[centre]
        symbol = symbols[centre]
        group = element(symbol).group_id
        square_planar = False
        if symbol in TRANSITION_METALS:
            n_lone_pairs = 0
            square_planar = (symbol, charges[centre]) in SQUARE_PLANAR
        else:
            valence = _valence_electrons(group)
            if valence is None or element(symbol).period < 3 and len(ligands) < 5:
                # Second-row main-group atoms are not hypervalent
                continue
            n_lone_pairs = max(
                0, int(round((valence - charges[centre] - order_sum[centre]) / 2))
            )
            if len(ligands) + n_lone_pairs < 5:
                continue
        if len(ligands) + n_lone_pairs not in _POLYHEDRA:
            continue

        branches = [_branch(neighbors, ligand, centre) for ligand in ligands]
        if any(
            other in branch
            for k, branch in enumerate(branches)
            for other in (ligands[:k] + ligands[k + 1 :])
        ):
            logger.warning(
                f"Atom {centre + 1} ({symbol}) has a chelating ligand; its geometry "
                "is left as the embedder made it."
            )
            continue

        targets = ligand_directions(len(ligands), n_lone_pairs, square_planar)
        # Match each ligand to the vertex nearest its present direction
        current = xyz[ligands] - xyz[centre]
        lengths = np.linalg.norm(current, axis=1)
        lengths[lengths < 1.0e-6] = 1.0
        current = current / lengths[:, None]
        free = list(range(len(targets)))
        assignment = {}
        for _, k, t in sorted(
            (
                (-float(current[k] @ targets[t]), k, t)
                for k in range(len(ligands))
                for t in free
            ),
        ):
            if k in assignment or t not in free:
                continue
            assignment[k] = t
            free.remove(t)

        new = xyz.copy()
        for k, ligand in enumerate(ligands):
            direction = targets[assignment[k]]
            donor = xyz[centre] + bond_length(symbol, symbols[ligand]) * direction
            rotation = _rotation(current[k], direction)
            for atom in branches[k]:
                new[atom] = donor + rotation @ (xyz[atom] - xyz[ligand])
        # Turn each ligand about its bond to the centre, to keep it clear of the
        # others (e.g. the methyls of P(OMe)6-)
        _turn_ligands(new, centre, ligands, branches)
        xyz = new
        changed.append(centre)
    return xyz, changed


def _axis_rotation(axis, angle):
    """The rotation matrix for ``angle`` (radians) about the unit vector ``axis``."""
    x, y, z = axis
    k = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * (k @ k)


def _turn_ligands(xyz, centre, ligands, branches, steps=24, passes=3):
    """Turn each polyatomic ligand about its bond to the centre to the angle that
    keeps its atoms farthest from those of the other ligands. Changes ``xyz``."""
    angles = np.linspace(0.0, 2.0 * np.pi, steps, endpoint=False)
    for _ in range(passes):
        for k, ligand in enumerate(ligands):
            moving = [a for a in branches[k] if a != ligand]
            if not moving:
                continue
            others = [centre] + [
                a for kk, branch in enumerate(branches) if kk != k for a in branch
            ]
            axis = xyz[ligand] - xyz[centre]
            axis /= np.linalg.norm(axis)
            relative = xyz[moving] - xyz[ligand]
            best_distance, best = -1.0, None
            for angle in angles:
                trial = xyz[ligand] + relative @ _axis_rotation(axis, angle).T
                d = np.linalg.norm(trial[:, None] - xyz[others][None], axis=2).min()
                if d > best_distance:
                    best_distance, best = d, trial
            xyz[moving] = best


def check_contacts(symbols, xyz, bonds, atoms=None):
    """Raise RuntimeError if atoms more than two bonds apart are too close, or any
    two fluorines are; only pairs involving ``atoms`` when given."""
    n = len(symbols)
    neighbors = [set() for _ in range(n)]
    for i, j, _ in bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
    xyz = np.asarray(xyz, dtype=float)
    selected = set(range(n)) if atoms is None else set(atoms)
    problems = []
    for i in range(n):
        for j in range(i + 1, n):
            if i not in selected or j not in selected:
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if j in neighbors[i]:
                continue
            if symbols[i] == "F" and symbols[j] == "F":
                if d < MIN_F_F:
                    problems.append((i, j, d, MIN_F_F))
                continue
            if neighbors[i] & neighbors[j]:
                continue  # 1-3 pairs are set by the bond angle
            limit = MIN_WITH_H if "H" in (symbols[i], symbols[j]) else MIN_HEAVY
            if d < limit:
                problems.append((i, j, d, limit))
    if problems:
        text = ", ".join(
            f"{symbols[i]}{i + 1}-{symbols[j]}{j + 1} {d:.2f} Å (< {limit})"
            for i, j, d, limit in problems[:5]
        )
        raise RuntimeError(f"The structure has atoms too close together: {text}.")

# -*- coding: utf-8 -*-

"""Topological methods for the system"""

import logging
from math import floor
import pprint  # noqa: F401

import numpy as np

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise

logger = logging.getLogger(__name__)

#: The charges of isolated atoms that are ions: alkali and alkaline-earth metals and
#: halogens, by atomic number
_ION_CHARGES = {
    **{z: 1 for z in (3, 11, 19, 37, 55)},
    **{z: 2 for z in (4, 12, 20, 38, 56)},
    **{z: -1 for z in (9, 17, 35, 53)},
}


#: The charges of the alkali and alkaline-earth metals in complexes, by atomic
#: number; the other metals share the rest of the charge
_METAL_ION_CHARGES = {
    **{z: 1 for z in (3, 11, 19, 37, 55, 87)},
    **{z: 2 for z in (4, 12, 20, 38, 56, 88)},
}


# Elements that are never bonded by distance-based perception unless asked:
# alkali and alkaline-earth metals, which are ions in molecular systems.
IONIC_ELEMENTS = ("Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba", "Ra")

# Default maximum number of bonds per element for perception; others unlimited.
MAX_BONDS = {"H": 1}

#: The metals for bond orders (assign_bond_orders): the alkali and alkaline-earth
#: metals, the transition metals, lanthanides and actinides, and the heavier
#: post-transition metals. B, Si, Ge, As, Sb, Se and Te are bonded covalently.
METALS = frozenset(
    (
        "Li Na K Rb Cs Fr Be Mg Ca Sr Ba Ra Al Ga In Tl Sn Pb Bi Po "
        "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
        "La Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs Mt Ds Rg Cn "
        "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
        "Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
    ).split()
)

#: The charges tried for a ligand, in order of preference among equals
_LIGAND_CHARGES = (0, -1, -2, 1, -3, 2)


def _ligand_lewis(atnos, xyzs, atoms, pairs, donors, all_options=False):
    """The Lewis structure of a ligand: a fragment of a metal complex without its
    metals, as RDKit's DetermineBondOrders gives it.

    Of the charges that give a structure without radicals, the one with the least
    charge separation (the smallest sum of |formal charges|) is taken, then the
    smallest |charge|, then a negative one: CO is neutral ([C-]#[O+]), Cp, Cl and
    OAc are -1, O is -2, H is -1. RDKit puts no lone pair on, e.g., the
    two-coordinate N of NMe2-, so a saturated ligand that has no structure
    otherwise gets its missing valence as negative charge on its donor atoms. A
    lone atom is the anion of its default valence: H-, Cl-, O2-, N3-, C4-.
    With `all_options`, all the structures found, as [(score, (orders, charges,
    charge))], best first.

    Parameters
    ----------
    atnos, xyzs : [int], [[float]]
        The atomic numbers and Cartesian coordinates of all the atoms.
    atoms : [int]
        The ligand's atoms (indices into atnos).
    pairs : [(int, int)]
        The ligand's bonds (indices into atnos).
    donors : set(int)
        Its atoms bonded to a metal.

    Returns
    -------
    ({(int, int): int}, {int: int}, int) or None
        The bond orders by pair, the formal charges by atom and the ligand's
        charge; or None if it has no closed-shell structure.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    to_order = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
        Chem.BondType.AROMATIC: 5,
    }
    table = Chem.GetPeriodicTable()
    if len(atoms) == 1:
        # A lone atom: the closed-shell anion of its default valence, e.g. H-,
        # Cl-, O2-, N3-, C4-, Si4-
        valence = max(0, table.GetDefaultValence(atnos[atoms[0]]))
        got = ({}, {atoms[0]: -valence}, -valence)
        return [((valence, valence, False), got)] if all_options else got
    local = {k: n for n, k in enumerate(atoms)}

    def molecule():
        mol = Chem.RWMol()
        for k in atoms:
            mol.AddAtom(Chem.Atom(atnos[k]))
        for i, j in pairs:
            mol.AddBond(local[i], local[j], Chem.BondType.SINGLE)
        conformer = Chem.Conformer(len(atoms))
        for k in atoms:
            conformer.SetAtomPosition(local[k], xyzs[k])
        mol.AddConformer(conformer, assignId=True)
        return mol.GetMol()

    def result(mol, q):
        orders = {}
        for i, j in pairs:
            bond = mol.GetBondBetweenAtoms(local[i], local[j])
            orders[(i, j)] = to_order.get(bond.GetBondType())
        if None in orders.values():
            return None
        charges = {k: mol.GetAtomWithIdx(local[k]).GetFormalCharge() for k in atoms}
        return orders, charges, q

    options = []
    for q in _LIGAND_CHARGES:
        mol = molecule()
        try:
            rdDetermineBonds.DetermineBondOrders(
                mol, charge=q, allowChargedFragments=True, embedChiral=False
            )
        except Exception:
            continue
        if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            continue
        score = (sum(abs(a.GetFormalCharge()) for a in mol.GetAtoms()), abs(q), q > 0)
        got = result(mol, q)
        if got is not None:
            options.append((score, got))
    if options:
        options.sort(key=lambda option: option[0])
        return options if all_options else options[0][1]

    # A saturated ligand whose donors lack valence: anions on the donors (NMe2-)
    mol = Chem.RWMol(molecule())
    q = 0
    for atom in mol.GetAtoms():
        missing = table.GetDefaultValence(atom.GetAtomicNum()) - atom.GetDegree()
        if missing > 0:
            if atoms[atom.GetIdx()] not in donors or missing > 2:
                return None
            atom.SetFormalCharge(-missing)
            q -= missing
        atom.SetNoImplicit(True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    got = result(mol, q)
    if all_options:
        return [((-q, abs(q), False), got)] if got is not None else []
    return got


def _balance_ligands(ligands, target, max_options=4, max_combinations=20000):
    """The ligands' structures, one each, whose charges add up to `target`, with
    the least charge separation in all; or None if there are none (or too many
    combinations to try).

    Parameters
    ----------
    ligands : [[(score, (orders, charges, q))]]
        Each ligand's structures, best first, as _ligand_lewis gives them.
    target : int
        The total charge of the ligands.
    """
    import itertools

    choices = [options[:max_options] for options in ligands]
    n = 1
    for options in choices:
        n *= len(options)
    if n > max_combinations:
        return None
    best = None
    for combination in itertools.product(*choices):
        if sum(got[2] for _, got in combination) != target:
            continue
        cost = tuple(sum(c) for c in zip(*(score for score, _ in combination)))
        if best is None or cost < best[0]:
            best = (cost, [got for _, got in combination])
    return None if best is None else best[1]


#: A metal's bond to the central atom of a four-membered chelate ring is dropped
#: when it is this much (Å) longer than its bonds to both chelating atoms
CHELATE_CENTER_MARGIN = 0.3


def _drop_chelate_centers(metal, Is, Js, Os, Ls):
    """Drop a metal's perceived bond to the central atom of a four-membered chelate
    ring: the C of a kappa-2 carboxylate, the N of a kappa-2 nitrate. The wider
    tolerance for metals (needed for pi ligands) reaches it, about 0.4 Å further
    than the two O, but it is no donor. A bond is dropped when the atom has two
    other neighbors bonded to the same metal, both by bonds at least
    CHELATE_CENTER_MARGIN shorter; a pi ligand's atoms are all at about the same
    distance from the metal (within 0.15 Å), so its bonds are kept.

    Parameters
    ----------
    metal : numpy.ndarray of bool
        Whether each atom is a metal.
    Is, Js, Os, Ls : [int], [int], [array], [float]
        The bonds' atoms, cell offsets and lengths (Å).

    Returns
    -------
    ([int], [int], [array])
        Is, Js and Os without the dropped bonds.
    """
    length = {}
    neighbors = {}
    for a, b, d in zip(Is, Js, Ls):
        length[(a, b)] = length[(b, a)] = d
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    drop = set()
    for k, (a, b, d) in enumerate(zip(Is, Js, Ls)):
        if metal[a] == metal[b]:
            continue
        m, x = (a, b) if metal[a] else (b, a)
        shorter = [
            y
            for y in neighbors[x]
            if y != m
            and not metal[y]
            and (m, y) in length
            and d >= length[(m, y)] + CHELATE_CENTER_MARGIN
        ]
        if len(shorter) >= 2:
            drop.add(k)
    if not drop:
        return Is, Js, Os
    keep = [k for k in range(len(Is)) if k not in drop]
    return [Is[k] for k in keep], [Js[k] for k in keep], [Os[k] for k in keep]


def _molecule_lewis(atnos, xyzs, fragment, pairs, charges, n_unpaired, _depth=0):
    """The Lewis structure of a molecule without metals from RDKit's
    DetermineBondOrders.

    For each of `charges`, a structure without formal charges is tried first and
    then one with them (N#As, not [N+]#[As-]). All hydrogens are explicit, so an
    atom short of its valence has radical electrons (RDKit would otherwise count
    the missing valence as implicit hydrogens). The first closed-shell structure
    is taken; failing that, the first whose radical electrons number
    `n_unpaired` (any number if None, none if 0). A doublet RDKit cannot build
    directly, such as the vinyl radical, is made from the closed-shell anion (or
    cation), its extra charge becoming the radical centre: C=[CH], not a
    triradical.

    Returns
    -------
    ({(int, int): int}, {int: int}, {int: int}, int) or None
        The bond orders by pair, the formal charges and radical electrons by atom,
        and the molecule's charge; or None.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    to_order = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
        Chem.BondType.AROMATIC: 5,
    }
    local = {k: n for n, k in enumerate(fragment)}

    def molecule():
        mol = Chem.RWMol()
        for k in fragment:
            mol.AddAtom(Chem.Atom(atnos[k]))
        for i, j in pairs:
            mol.AddBond(local[i], local[j], Chem.BondType.SINGLE)
        conformer = Chem.Conformer(len(fragment))
        for k in fragment:
            conformer.SetAtomPosition(local[k], xyzs[k])
        mol.AddConformer(conformer, assignId=True)
        return mol.GetMol()

    second = None
    for q in charges:
        for charged in (False, True):
            mol = molecule()
            try:
                rdDetermineBonds.DetermineBondOrders(
                    mol, charge=q, allowChargedFragments=charged
                )
            except Exception:
                continue
            if Chem.GetFormalCharge(mol) != q:
                continue
            orders = {}
            for i, j in pairs:
                bond = mol.GetBondBetweenAtoms(local[i], local[j])
                orders[(i, j)] = to_order.get(bond.GetBondType())
            if None in orders.values():
                continue
            # The radicals with every hydrogen explicit
            for atom in mol.GetAtoms():
                atom.SetNoImplicit(True)
                atom.SetNumRadicalElectrons(0)
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                continue
            charges_ = {
                k: mol.GetAtomWithIdx(local[k]).GetFormalCharge() for k in fragment
            }
            radicals = {
                k: mol.GetAtomWithIdx(local[k]).GetNumRadicalElectrons()
                for k in fragment
            }
            n_radicals = sum(radicals.values())
            if n_radicals == 0:
                return orders, charges_, radicals, q
            if second is None and n_unpaired != 0:
                if n_unpaired is None or n_radicals == n_unpaired:
                    second = (orders, charges_, radicals, q)
    if second is not None or n_unpaired != 1 or _depth > 0:
        return second

    # A doublet RDKit cannot build directly (the vinyl radical): the closed-shell
    # anion (or cation), with its extra charge made the radical centre
    for shift, pick in ((-1, min), (1, max)):
        ion = _molecule_lewis(
            atnos,
            xyzs,
            fragment,
            pairs,
            [q + shift for q in charges],
            0,
            _depth=1,
        )
        if ion is None:
            continue
        orders, charges_, radicals, q = ion
        centers = [k for k in fragment if charges_[k] * shift > 0]
        # The least electronegative anionic (most electronegative cationic) atom
        # first, if it can hold the radical: a neutral B with four bonds cannot
        centers.sort(key=lambda k: _ELECTRONEGATIVITY.get(atnos[k], 2.5))
        if pick is max:
            centers.reverse()
        for k in centers:
            trial_charges = dict(charges_)
            trial_radicals = dict(radicals)
            trial_charges[k] -= shift
            trial_radicals[k] += 1
            if _valid(atnos, fragment, orders, trial_charges, trial_radicals):
                return orders, trial_charges, trial_radicals, q - shift
    return None


def _valid(atnos, fragment, orders, charges, radicals):
    """Whether RDKit accepts the structure, with every hydrogen explicit."""
    from rdkit import Chem

    to_type = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        5: Chem.BondType.AROMATIC,
    }
    local = {k: n for n, k in enumerate(fragment)}
    mol = Chem.RWMol()
    for k in fragment:
        atom = Chem.Atom(atnos[k])
        atom.SetFormalCharge(charges[k])
        atom.SetNumRadicalElectrons(radicals[k])
        atom.SetNoImplicit(True)
        mol.AddAtom(atom)
    for (i, j), order in orders.items():
        mol.AddBond(local[i], local[j], to_type[order])
        if order == 5:
            bond = mol.GetBondBetweenAtoms(local[i], local[j])
            bond.SetIsAromatic(True)
    mol = mol.GetMol()
    try:
        mol.UpdatePropertyCache(strict=True)
    except Exception:
        return False
    # The radicals RDKit would assign must be the ones given
    for k in fragment:
        atom = mol.GetAtomWithIdx(local[k])
        if atom.GetNumRadicalElectrons() != radicals[k]:
            return False
    return True


#: Pauling electronegativities, by atomic number, to choose the radical centre
_ELECTRONEGATIVITY = {
    1: 2.20,
    5: 2.04,
    6: 2.55,
    7: 3.04,
    8: 3.44,
    9: 3.98,
    14: 1.90,
    15: 2.19,
    16: 2.58,
    17: 3.16,
    32: 2.01,
    33: 2.18,
    34: 2.55,
    35: 2.96,
    50: 1.96,
    51: 2.05,
    52: 2.10,
    53: 2.66,
}


def _radicals_of_single_bonds(atnos, fragment, pairs):
    """The radical electrons of a neutral molecule with only single bonds, as RDKit
    assigns them (each atom's lowest allowed valence), or None if RDKit does not
    accept it (e.g. hypervalent PF6)."""
    from rdkit import Chem

    local = {k: n for n, k in enumerate(fragment)}
    mol = Chem.RWMol()
    for k in fragment:
        atom = Chem.Atom(atnos[k])
        atom.SetNoImplicit(True)
        mol.AddAtom(atom)
    for i, j in pairs:
        mol.AddBond(local[i], local[j], Chem.BondType.SINGLE)
    mol = mol.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return {k: mol.GetAtomWithIdx(local[k]).GetNumRadicalElectrons() for k in fragment}


def _components(atoms, pairs):
    """The connected components of the graph of `atoms` joined by `pairs`."""
    neighbors = {k: [] for k in atoms}
    for i, j in pairs:
        neighbors[i].append(j)
        neighbors[j].append(i)
    seen = set()
    result = []
    for start in atoms:
        if start in seen:
            continue
        seen.add(start)
        component = [start]
        stack = [start]
        while stack:
            k = stack.pop()
            for n in neighbors[k]:
                if n not in seen:
                    seen.add(n)
                    component.append(n)
                    stack.append(n)
        result.append(sorted(component))
    return result


class TopologyMixin:
    """A mixin for handling topology in a configuration."""

    def perceive_bonds(
        self,
        method="covalent radii",
        replace=False,
        tolerance=1.2,
        metal_tolerance=1.3,
        radii=None,
        max_bonds=None,
        exclude=IONIC_ELEMENTS,
        min_distance=0.4,
    ):
        """Perceive the bonds from the geometry and add them to the configuration.

        Bonds are single bonds; no attempt is made to assign bond orders. For a
        periodic configuration every periodic image within reach is considered
        and each bond is stored with the cell offset of its second atom (the
        CIF-style ``symop2`` such as ``1_556``), so molecules that straddle the
        cell boundary are bonded correctly and covalent crystals get all their
        bonds -- including several bonds between the same two atoms through
        different images (primitive diamond) and bonds from an atom to its own
        image (a one-atom cell). Any cell shape is handled. Symmetry other than
        P1 is not supported.

        Parameters
        ----------
        method : str = "covalent radii"
            The perception method. Currently only "covalent radii": atoms i and j
            are bonded if their distance is less than ``tolerance`` times the sum
            of their covalent radii. Other methods (e.g. Voronoi) may be added.
        replace : bool = False
            Whether to delete any existing bonds first. If False and the
            configuration already has bonds, an error is raised.
        tolerance : float = 1.2
            The factor multiplying the sum of the covalent radii.
        metal_tolerance : float = 1.3
            The factor for pairs with a metal (:data:`METALS`), whose bonds to
            pi ligands (Cp, arenes, dienes) are longer than covalent radii suggest:
            at 1.2 the Cr-C bonds of (C6H3Me3)Cr(CO)3 are missed. None uses
            ``tolerance``.
        radii : {str: float} = None
            Covalent radii in Å, by element symbol, overriding the defaults
            (Pyykkö radii from ``mendeleev``) for those elements.
        max_bonds : {str: int} = None
            The maximum number of bonds per atom, by element symbol, overriding
            the defaults (H: 1). When an atom would exceed its limit the shortest
            candidate bonds are kept.
        exclude : iterable of str = IONIC_ELEMENTS
            Elements that never form bonds -- by default the alkali and
            alkaline-earth metals, which are ions in molecular systems. Pass an
            empty tuple to bond them (e.g. for ionic crystals).
        min_distance : float = 0.4
            Pairs closer than this (Å) are ignored rather than bonded, guarding
            against duplicate or overlapping atoms.

        Returns
        -------
        int
            The number of bonds added.
        """
        if method != "covalent radii":
            raise ValueError(f"Unknown bond perception method '{method}'")
        if self.symmetry.n_symops > 1:
            raise NotImplementedError(
                "Bond perception is only supported for configurations without "
                "symmetry (P1 / C1)."
            )
        bonds = self.bonds
        if bonds.n_bonds > 0:
            if replace:
                bonds.delete()
            else:
                raise RuntimeError(
                    "The configuration already has bonds. Use replace=True to "
                    "perceive them afresh."
                )

        n_atoms = self.atoms.n_atoms
        if n_atoms == 0 or (n_atoms < 2 and self.periodicity == 0):
            return 0

        symbols = self.atoms.symbols
        unique = sorted(set(symbols))
        radius = self._covalent_radii(unique)
        if radii is not None:
            radius.update({k: float(v) for k, v in radii.items()})
        limits = dict(MAX_BONDS)
        if max_bonds is not None:
            limits.update(max_bonds)
        excluded = set(exclude) if exclude is not None else set()

        r_atom = np.array([radius[s] for s in symbols])
        active = np.array([s not in excluded for s in symbols])
        metal = np.array([s in METALS for s in symbols])
        if metal_tolerance is None:
            metal_tolerance = tolerance
        widest = max(tolerance, metal_tolerance) if metal.any() else tolerance
        cutoff = widest * 2.0 * r_atom[active].max() if active.any() else 0.0
        if cutoff <= 0.0:
            return 0

        pairs, dist, offsets = self._neighbor_pairs(cutoff)

        # Apply the criterion pair by pair.
        i, j = pairs[:, 0], pairs[:, 1]
        keep = (
            active[i]
            & active[j]
            & (dist > min_distance)
            & (
                dist
                < np.where(metal[i] | metal[j], metal_tolerance, tolerance)
                * (r_atom[i] + r_atom[j])
            )
        )
        i, j, dist, offsets = i[keep], j[keep], dist[keep], offsets[keep]

        # Enforce the per-element bond limits, keeping the shortest bonds.
        order = np.argsort(dist, kind="stable")
        count = np.zeros(n_atoms, dtype=int)
        unlimited = 10**9
        limit = np.array([limits.get(s, unlimited) for s in symbols])
        Is, Js, Os, Ls = [], [], [], []
        for k in order:
            a, b = int(i[k]), int(j[k])
            if a == b:
                # A bond from an atom to its own periodic image counts twice.
                if count[a] + 2 > limit[a]:
                    continue
                count[a] += 2
            else:
                if count[a] >= limit[a] or count[b] >= limit[b]:
                    continue
                count[a] += 1
                count[b] += 1
            Is.append(a)
            Js.append(b)
            Os.append(offsets[k])
            Ls.append(float(dist[k]))

        if self.periodicity == 0 and metal.any():
            Is, Js, Os = _drop_chelate_centers(metal, Is, Js, Os, Ls)

        if len(Is) == 0:
            return 0

        ids = self.atoms.ids
        kwargs = {}
        if self.periodicity != 0:
            # Cell offsets of the second atom, in the CIF-style '1_555' encoding
            # (identity operator; each digit is 5 + the offset along that axis).
            if np.abs(np.array(Os)).max() > 4:
                raise ValueError(
                    "A perceived bond spans more than 4 cell repeats; the cell is "
                    "too small to encode the bond offsets."
                )
            kwargs["symop1"] = ["."] * len(Is)
            kwargs["symop2"] = [
                "." if not np.any(o) else "1_" + "".join(str(5 + int(x)) for x in o)
                for o in Os
            ]
        bonds.append(
            i=[ids[a] for a in Is],
            j=[ids[b] for b in Js],
            bondorder=[1] * len(Is),
            **kwargs,
        )
        return len(Is)

    def assign_bond_orders(self, charge=None):
        """Assign bond orders to the perceived (single) bonds, molecule by molecule.

        :meth:`perceive_bonds` gives connectivity only, so ethylene would be read
        as [CH2+][CH2-] by anything that derives a SMILES or a Lewis structure.
        Each molecule (connected fragment) without metals gets the bond orders,
        formal charges and radical electrons of RDKit's ``DetermineBondOrders``:
        the first closed-shell structure of the charges tried (the molecule's own
        charge when it is the only one, else 0, -1, +1, -2, +2), or failing that
        one with as many radical electrons as the spin multiplicity says (the
        ethyl or vinyl radical). A molecule with neither keeps its single bonds and
        gets radical electrons where its atoms lack valence (Ge2 as [Ge][Ge], GeC
        as [C][Ge], not charge-separated), when that keeps the net charge;
        otherwise (hypervalent PF6-, whose single bonds are already its Lewis
        structure) it keeps its bonds as they are. Single atoms of the alkali and
        alkaline-earth metals and of the halogens get their usual ionic charges.

        Metal complexes (molecules with an atom in :data:`METALS`) are ionic +
        dative: each ligand, a fragment of the molecule without its metals, gets
        its closed-shell Lewis structure with the least charge separation (CO 0,
        Cp-, Cl-, OAc-, NMe2-, O2-, H-, C4-; see :func:`_ligand_lewis`); every
        metal-ligand bond is a dative bond from the ligand's donor atom to the
        metal (:data:`molsystem.bonds.DATIVE`); metal-metal bonds are single
        bonds; and the metals carry the rest of the net charge: alkali metals +1
        and alkaline-earth metals +2, the others sharing what is left as evenly
        as integers allow. The unpaired electrons of an open-shell complex are
        the metals'. The metals' formal charges are bookkeeping, not oxidation
        states. Long metal-ligand contacts, agostic hydrogens, bridging hydrides
        and clusters depend on which bonds were perceived and are not checked.

        Periodic configurations are left alone: a cell's bonds cross the cell.

        Parameters
        ----------
        charge : int = None
            The net charge, if known, to check the molecules' charges against. The
            metals of complexes take what is left of it (of the configuration's
            charge if not given).

        Returns
        -------
        (int, int)
            The number of molecules assigned and the number that could not be. If
            the molecules' charges do not add up to `charge`, a warning is logged.
        """
        from rdkit import Chem, RDLogger

        from .bonds import DATIVE

        if self.periodicity != 0 or self.bonds.n_bonds == 0:
            return 0, 0

        ids = list(self.atoms.ids)
        index = {_id: k for k, _id in enumerate(ids)}
        rows = [(index[row["i"]], index[row["j"]]) for row in self.bonds.bonds()]
        atnos = [int(a) for a in self.atoms.atomic_numbers]
        symbols = list(self.atoms.symbols)
        is_metal = [symbol in METALS for symbol in symbols]
        xyzs = [
            [float(x) for x in xyz]
            for xyz in self.atoms.get_coordinates(fractionals=False)
        ]

        whole = Chem.RWMol()
        for atno in atnos:
            whole.AddAtom(Chem.Atom(atno))
        for i, j in rows:
            whole.AddBond(i, j, Chem.BondType.SINGLE)
        fragments = Chem.GetMolFrags(whole.GetMol(), sanitizeFrags=False)
        complexes = [f for f in fragments if len(f) > 1 and any(is_metal[k] for k in f)]

        if charge is None:
            charge = getattr(self, "charge", None)
        multiplicity = getattr(self, "spin_multiplicity", None)
        n_unpaired = None if multiplicity is None else multiplicity - 1
        only_one = len(fragments) == 1

        order = {pair: 1 for pair in rows}
        formal = [0] * len(atnos)
        if "formal_charge" in self.atoms:
            formal = [int(q or 0) for q in self.atoms.get_column_data("formal_charge")]
        radicals = [0] * len(atnos)
        n_ok = n_failed = 0
        total = 0
        failed = []  # molecules with no Lewis structure, for radicals afterwards
        assigned = []  # molecules without metals given a structure
        RDLogger.DisableLog("rdApp.*")
        try:
            for fragment in fragments:
                if len(fragment) == 1:
                    k = fragment[0]
                    ion = _ION_CHARGES.get(atnos[k])
                    if ion is not None:
                        formal[k] = ion
                        total += ion
                        n_ok += 1
                    else:
                        total += formal[k]
                    continue
                if fragment in complexes:
                    continue  # after the other molecules, for the metals' charges
                pairs = [(i, j) for i, j in rows if i in set(fragment)]
                if only_one and charge is not None:
                    charges = (charge,)
                else:
                    charges = (0, -1, 1, -2, 2)
                if only_one:
                    expected = n_unpaired
                else:
                    expected = 0 if n_unpaired == 0 else None
                result = _molecule_lewis(
                    atnos, xyzs, fragment, pairs, charges, expected
                )
                if result is None:
                    failed.append((fragment, pairs))
                    continue
                orders, charges_, radicals_, q = result
                order.update(orders)
                for k in fragment:
                    formal[k] = charges_[k]
                    radicals[k] = radicals_[k]
                total += q
                n_ok += 1
                assigned.append((fragment, pairs, q))

            # Ions and molecules whose charges do not add up (Cs2O2: Cs+ x 2 and
            # O=O): other closed-shell charges of the molecules that balance them
            # (peroxide); failing that, the metal ions give way (SrOH: Sr+, OH-)
            # Only when every molecule has a structure: an unassigned one (PF6-)
            # carries the difference
            if not complexes and not failed and charge is not None and total != charge:
                total = self._balance_ions(
                    fragments,
                    assigned,
                    atnos,
                    xyzs,
                    is_metal,
                    charge,
                    total,
                    order,
                    formal,
                    radicals,
                )

            # The metal complexes: ligands, dative bonds, and the metals' charges
            metals = []
            ligands = []  # each ligand's closed-shell structures, best first
            for fragment in complexes:
                in_fragment = set(fragment)
                metals.extend(k for k in fragment if is_metal[k])
                ligand_pairs = [
                    (i, j)
                    for i, j in rows
                    if i in in_fragment and not is_metal[i] and not is_metal[j]
                ]
                donors = set()
                for i, j in rows:
                    if i in in_fragment and is_metal[i] != is_metal[j]:
                        donors.add(j if is_metal[i] else i)
                        order[(i, j)] = DATIVE if is_metal[j] else -DATIVE
                    elif i in in_fragment and is_metal[i] and is_metal[j]:
                        order[(i, j)] = 1
                ligand_failed = False
                for ligand in _components(
                    [k for k in fragment if not is_metal[k]], ligand_pairs
                ):
                    pairs = [(i, j) for i, j in ligand_pairs if i in ligand]
                    options = _ligand_lewis(
                        atnos, xyzs, ligand, pairs, donors, all_options=True
                    )
                    if not options:
                        ligand_failed = True
                        total += sum(formal[k] for k in ligand)
                        continue
                    ligands.append(options)
                if ligand_failed:
                    n_failed += 1
                    logger.warning(
                        "A ligand of a metal complex has no closed-shell Lewis "
                        "structure; its bonds are left single."
                    )
                else:
                    n_ok += 1

            # Molecules with no Lewis structure: single bonds and radicals, if
            # that keeps the net charge (it does not for PF6-)
            if failed:
                rest = None if charge is None or metals else charge - total
                for fragment, pairs in failed:
                    found = None
                    if rest == 0 or (rest is None and charge is None):
                        found = _radicals_of_single_bonds(atnos, fragment, pairs)
                    if found is None:
                        n_failed += 1
                        total += sum(formal[k] for k in fragment)
                        continue
                    for k in fragment:
                        formal[k] = 0
                        radicals[k] = found[k]
                    n_ok += 1

            if metals:
                fixed = {k: _METAL_ION_CHARGES.get(atnos[k]) for k in metals}
                others = [k for k in metals if fixed[k] is None]
                # Each ligand's best structure; with only alkali and
                # alkaline-earth metals, whose charges are fixed, the ones that
                # balance the net charge if there are (Cs2O2: peroxide, not O2)
                chosen = [options[0][1] for options in ligands]
                if not others and charge is not None and ligands:
                    target = charge - total - sum(fixed.values())
                    if sum(got[2] for got in chosen) != target:
                        balanced = _balance_ligands(ligands, target)
                        if balanced is not None:
                            chosen = balanced
                for orders, charges_, q in chosen:
                    order.update(orders)
                    for k, value in charges_.items():
                        formal[k] = value
                        radicals[k] = 0
                    total += q
                rest = (charge or 0) - total
                if others:
                    for k, q in fixed.items():
                        if q is not None:
                            formal[k] = q
                            rest -= q
                    takers = others
                elif charge is not None and rest != sum(fixed.values()):
                    # The fixed charges give way when the ligands cannot balance
                    # them (SrOH: Sr+)
                    takers = metals
                else:
                    takers = []
                    for k, q in fixed.items():
                        formal[k] = q
                    rest = 0
                if takers:
                    share = [rest // len(takers)] * len(takers)
                    for k in range(rest - sum(share)):
                        share[k] += 1
                    for k, q in zip(takers, share):
                        formal[k] = q
                total = sum(formal[k] for k in range(len(atnos)))
                for k in metals:
                    radicals[k] = 0
        finally:
            RDLogger.EnableLog("rdApp.*")

        if n_ok > 0:
            Is, Js, Os = [], [], []
            for i, j in rows:
                bondorder = order[(i, j)]
                if bondorder == -DATIVE:  # the donor is j
                    i, j, bondorder = j, i, DATIVE
                Is.append(ids[i])
                Js.append(ids[j])
                Os.append(bondorder)
            self.bonds.delete()
            self.bonds.append(i=Is, j=Js, bondorder=Os)
            if "formal_charge" not in self.atoms:
                self.atoms.add_attribute("formal_charge", coltype="int", default=0)
            self.atoms["formal_charge"][0:] = formal
            # Recorded even when all zero: the structure assigned here is the one
            # to keep (to_RDKMol uses it rather than fitting the multiplicity)
            if "radical_electrons" not in self.atoms:
                self.atoms.add_attribute("radical_electrons", coltype="int", default=0)
            self.atoms["radical_electrons"][0:] = radicals
        if charge is not None and total != charge:
            logger.warning(
                f"The molecules' charges add up to {total}, not the net charge "
                f"{charge}: check the bond orders."
            )
        return n_ok, n_failed

    @staticmethod
    def _balance_ions(
        fragments,
        assigned,
        atnos,
        xyzs,
        is_metal,
        charge,
        total,
        order,
        formal,
        radicals,
    ):
        """Make the charges of ions and molecules add up to the net charge: other
        closed-shell charges of the molecules, with the least charge separation,
        if there are some that balance (Cs2O2: O2 2-); otherwise the single metal
        ions take what is left (SrOH: Sr+). Updates order, formal and radicals in
        place and returns the new total."""
        ions = [f[0] for f in fragments if len(f) == 1 and is_metal[f[0]]]
        if assigned:
            others = total - sum(q for _, _, q in assigned)
            options = []
            for fragment, pairs, _ in assigned:
                found = []
                for rank, q in enumerate((0, -1, 1, -2, 2)):
                    got = _molecule_lewis(atnos, xyzs, fragment, pairs, (q,), 0)
                    if got is not None:
                        found.append(((rank,), got))
                options.append(found)
            if all(options):
                chosen = _balance_ligands(
                    [
                        [
                            (score, (got[0], got[1], got[3], got[2]))
                            for score, got in found
                        ]
                        for found in options
                    ],
                    charge - others,
                )
                if chosen is not None:
                    for (fragment, _, _), (orders, charges_, q, radicals_) in zip(
                        assigned, chosen
                    ):
                        order.update(orders)
                        for k in fragment:
                            formal[k] = charges_[k]
                            radicals[k] = radicals_[k]
                    total = sum(formal)
                    if total == charge:
                        return total
        if ions:
            rest = charge - (total - sum(formal[k] for k in ions))
            share = [rest // len(ions)] * len(ions)
            for k in range(rest - sum(share)):
                share[k] += 1
            for k, q in zip(ions, share):
                formal[k] = q
            total = sum(formal)
        return total

    @staticmethod
    def _covalent_radii(symbols):
        """Pyykkö covalent radii (Å) for the given element symbols, via mendeleev."""
        import mendeleev

        result = {}
        for symbol in symbols:
            element = mendeleev.element(symbol)
            r = element.covalent_radius_pyykko
            if r is None:
                r = element.covalent_radius
            if r is None:
                raise ValueError(f"No covalent radius available for {symbol}")
            result[symbol] = r / 100.0  # pm -> Å
        return result

    def _neighbor_pairs(self, cutoff):
        """All atom pairs within ``cutoff`` Å, each periodic image separately.

        For a periodic configuration every image of atom j within the cutoff of
        atom i is a distinct candidate bond, so a small cell (e.g. primitive
        diamond) gets all its bonds, including bonds from an atom to its own
        image. Each bond is reported once: pairs have i <= j, and a self-image
        bond is reported for one of the two opposite offsets only.

        Returns
        -------
        (pairs, distances, offsets)
            pairs : ndarray (n, 2) of 0-based atom indices, i <= j
            distances : ndarray (n,) in Å
            offsets : ndarray (n, 3) of int -- the cell offset of atom j
                relative to the stored coordinates (all zero if not periodic)
        """
        from scipy.spatial import cKDTree

        empty = (np.zeros((0, 2), dtype=int), np.zeros(0), np.zeros((0, 3), dtype=int))
        periodicity = self.periodicity
        if periodicity == 0:
            xyz = self.atoms.get_coordinates(fractionals=False, as_array=True)
            pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
            if len(pairs) == 0:
                return empty
            pairs = np.sort(pairs, axis=1)
            dist = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
            return pairs, dist, np.zeros((len(pairs), 3), dtype=int)

        if periodicity != 3:
            raise NotImplementedError(
                "Bond perception is only implemented for 3-D periodic and "
                "non-periodic configurations."
            )

        T = self.cell.to_cartesians_transform(as_array=True)  # xyz = uvw @ T
        stored = self.atoms.get_coordinates(fractionals=True, as_array=True)
        wrap = np.floor(stored).astype(int)
        uvw = stored - wrap  # in [0, 1)
        n_atoms = len(uvw)

        # How many image shells are needed along each axis: the cutoff divided
        # by the spacing between the lattice planes perpendicular to that axis.
        volume = abs(np.linalg.det(T))
        spacing = [
            volume / np.linalg.norm(np.cross(T[(k + 1) % 3], T[(k + 2) % 3]))
            for k in range(3)
        ]
        shells = [int(np.ceil(cutoff / d)) for d in spacing]
        ranges = [range(-n, n + 1) for n in shells]
        shifts = np.array(
            [[a, b, c] for a in ranges[0] for b in ranges[1] for c in ranges[2]]
        )

        central = uvw @ T
        images = (uvw[None, :, :] + shifts[:, None, :]).reshape(-1, 3) @ T
        image_atom = np.tile(np.arange(n_atoms), len(shifts))
        image_shift = np.repeat(shifts, n_atoms, axis=0)

        tree = cKDTree(images)
        hits = tree.query_ball_point(central, cutoff)
        Is, Js, D, Offs = [], [], [], []
        for i, neigh in enumerate(hits):
            for h in neigh:
                j = int(image_atom[h])
                if j < i:
                    continue  # reported from j's side with the opposite shift
                s = image_shift[h]
                if j == i:
                    if not np.any(s):
                        continue  # the atom itself
                    # keep one of the two opposite offsets
                    if tuple(s) < tuple(-s):
                        continue
                r = float(np.linalg.norm(images[h] - central[i]))
                Is.append(i)
                Js.append(j)
                D.append(r)
                # offset of j relative to the *stored* coordinates of i and j
                Offs.append(s - wrap[j] + wrap[i])
        if len(Is) == 0:
            return empty
        pairs = np.stack([np.array(Is), np.array(Js)], axis=1)
        return pairs, np.array(D), np.array(Offs, dtype=int)

    def find_molecules(self, as_indices=False):
        """Find the separate molecules.

        Parameters
        ----------
        as_indices : bool = False
            Whether to return 0-based indices (True) or atom ids (False)

        Returns
        -------
        molecules : [[int]*n_molecules]
            A list of lists of atom ids or indices for the molecules
        """

        molecules = []

        n_atoms = self.atoms.n_atoms

        if n_atoms == 0:
            return molecules

        if not as_indices and self.symmetry.n_symops > 1:
            raise RuntimeError(
                "Cannot return atom ids for bonded_neighbors when there is symmetry"
            )
        neighbors = self.bonded_neighbors(as_indices=True)
        visited = [False] * n_atoms
        while True:
            # Find first atom not yet visited
            try:
                i = visited.index(False)
            except ValueError:
                break
            visited[i] = True
            atoms = [i]
            next_atoms = neighbors[i]
            while len(next_atoms) > 0:
                tmp = []
                for i in next_atoms:
                    if not visited[i]:
                        atoms.append(i)
                        visited[i] = True
                        tmp.extend(neighbors[i])
                next_atoms = tmp
            molecules.append(sorted(atoms))
        if as_indices:
            return molecules
        else:
            atom_ids = self.atoms.ids
            to_id = {i: j for i, j in enumerate(atom_ids)}
            return [[to_id[j] for j in js] for js in molecules]

    def bonded_neighbors(self, as_indices=False, first_index=0):
        """The atoms bonded to each atom in the system.

        Parameters
        ----------
        as_indices : bool = False
            Whether to return 0-based indices (True) or atom ids (False)
        first_index : int = 0
            The smallest index, e.g. 0 or 1

        Returns
        -------
        neighbors : {int: [int]} or [[int]] for indices
            list of atom ids for each atom id
        """
        neighbors = {}

        symmetry = self.symmetry
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = atoms.n_atoms

        if n_atoms == 0:
            if as_indices:
                return []
            else:
                return neighbors

        if self.symmetry.n_symops > 1:
            # Symmetry involved...
            if not as_indices:
                raise RuntimeError(
                    "Cannot return atom ids for bonded_neighbors when there is symmetry"
                )
            pairs = symmetry.bond_atoms
            neighbors = [[] for i in range(n_atoms)]
            for i, j in pairs:
                neighbors[i].append(j)
                neighbors[j].append(i)
            for js in neighbors:
                js.sort()
            return neighbors
        else:
            atom_ids = atoms.ids
            neighbors = {i: [] for i in atom_ids}

            if bonds.n_bonds > 0:
                for bond in bonds.bonds():
                    i = bond["i"]
                    j = bond["j"]
                    neighbors[i].append(j)
                    neighbors[j].append(i)

            if as_indices:
                # Convert to indices
                to_index = {j: i + first_index for i, j in enumerate(atom_ids)}
                result = [[]] * (n_atoms + first_index)
                for i, js in neighbors.items():
                    result[to_index[i]] = sorted([to_index[j] for j in js])
                return result
            else:
                for i in neighbors:
                    neighbors[i].sort()

                return neighbors

    def create_molecule_subsets(self):
        """Create a subset for each molecule in a configuration.

        Returns
        -------
        [int]
            The ids of the subsets, one per molecule.
        """
        # Find the molecules and the create the subsets if they don't exist.
        molecules = self.find_molecules()

        # Remove any previous subsets for this configuration
        subsets = self["subset"]
        tid = 1
        sids = subsets.find(tid)
        if len(sids) > 0:
            subsets.delete(sids)

        # Now create the new set.
        sids = []
        for atom_ids in molecules:
            sid = subsets.create(tid, atoms=atom_ids)
            sids.append(sid)

        return sids

    def create_molecule_templates(self, full_templates=True, create_subsets=True):
        """Create a template for each unique molecule in a configuration.

        By default also create subsets linking each template to the atoms
        of the molecules in the system.

        Parameters
        ----------
        full_templates : bool = True
            If true, create full templates by creating systems for the
            molecules.
        create_subsets : bool = True
            If true, create subsets linking the templates to the molecules.

        Returns
        -------
        [int] or [[int], [int]]
            The ids of the templates, or if create_subsets is True
            a two-element list containing the list of templates and
            list of subsets.
        """
        templates = self.system_db.templates

        # Find the molecules
        molecules = self.find_molecules()
        n_molecules = len(molecules)

        # And the molecule each atom is in
        atom_to_molecule = {}
        for molecule, atoms in enumerate(molecules):
            for atom in atoms:
                atom_to_molecule[atom] = molecule

        # The bonds in each molecule
        bonds_per_molecule = [[] for i in range(n_molecules)]
        for bond in self.bonds.bonds():
            i = bond["i"]
            j = bond["j"]
            order = bond["bondorder"]
            molecule = atom_to_molecule[i]
            bonds_per_molecule[molecule].append((i, j, order))

        # Get the canonical smiles for each molecule
        to_can = openbabel.OBConversion()
        to_can.SetOutFormat("can")
        ob_mol = openbabel.OBMol()
        ob_template = openbabel.OBMol()
        atnos = self.atoms.atomic_numbers
        xyzs = self.atoms.get_coordinates(fractionals=False)

        atom_index = {j: i for i, j in enumerate(self.atoms.ids)}

        new_subsets = {}
        sids = {}
        new_templates = []
        tids = []
        for molecule, atoms in enumerate(molecules):
            to_index = {j: i for i, j in enumerate(atoms)}
            molecule_atnos = [atnos[atom_index[i]] for i in atoms]

            ob_mol.Clear()
            for atom, atno in zip(atoms, molecule_atnos):
                ob_atom = ob_mol.NewAtom()
                ob_atom.SetAtomicNum(atno)
            bonds = []
            for i, j, order in bonds_per_molecule[molecule]:
                bonds.append((to_index[i], to_index[j], order))
                # 1-based indices in ob.
                ob_mol.AddBond(to_index[i] + 1, to_index[j] + 1, order)

            canonical = to_can.WriteString(ob_mol).strip()

            if full_templates:
                # See if a molecule template with the canonical smiles exists
                if templates.exists(canonical, "molecule"):
                    template = templates.get(canonical, category="molecule")
                else:
                    # Create a new system & configuration for the template
                    system_name = "template system " + canonical
                    if not self.system_db.system_exists(system_name):
                        system = self.system_db.create_system(
                            system_name, make_current=False
                        )
                        configuration = system.create_configuration(
                            canonical, make_current=False
                        )
                        cid = configuration.id

                        kwargs = {}
                        kwargs["atno"] = molecule_atnos
                        molecule_xyzs = [xyzs[atom_index[i]] for i in atoms]
                        kwargs["x"] = [x for x, y, z in molecule_xyzs]
                        kwargs["y"] = [y for x, y, z in molecule_xyzs]
                        kwargs["z"] = [z for x, y, z in molecule_xyzs]

                        ids = configuration.atoms.append(**kwargs)

                        kwargs = {}
                        kwargs["i"] = [ids[x] for x, _, _ in bonds]
                        kwargs["j"] = [ids[x] for _, x, _ in bonds]
                        kwargs["bondorder"] = [x for _, _, x in bonds]

                        configuration.bonds.append(**kwargs)

                        template = templates.create(
                            canonical, category="molecule", configuration=cid
                        )
            else:
                if templates.exists(canonical, "molecule"):
                    template = templates.get(canonical, category="molecule")
                else:
                    template = templates.create(canonical, category="molecule")

            if template.id not in tids:
                tids.append(template.id)
                new_templates.append(template)

            if create_subsets:
                if full_templates:
                    # Need to reorder the atoms to match the template atoms
                    # Prepare the OB molecule for the template
                    ob_template.Clear()
                    for atno in template.atoms.atomic_numbers:
                        ob_atom = ob_template.NewAtom()
                        ob_atom.SetAtomicNum(atno)

                    tatom_ids = template.atoms.ids
                    to_index = {j: i for i, j in enumerate(tatom_ids)}
                    for row in template.bonds.bonds():
                        i = to_index[row["i"]]
                        j = to_index[row["j"]]
                        order = row["bondorder"]
                        ob_template.AddBond(i + 1, j + 1, order)

                    # Get the mapping from template to molecule
                    query = openbabel.CompileMoleculeQuery(ob_template)
                    mapper = openbabel.OBIsomorphismMapper.GetInstance(query)
                    mapping = openbabel.vpairUIntUInt()
                    mapper.MapFirst(ob_mol, mapping)
                    reordered_atoms = [atoms[j] for i, j in mapping]

                    if len(reordered_atoms) == 0:
                        logger.warning(
                            "There are no atoms in the subset for this molecule!"
                        )
                    else:
                        subset = self.subsets.create(
                            template, reordered_atoms, tatom_ids
                        )
                else:
                    subset = self.subsets.create(template, atoms)

                tid = template.id
                if tid not in sids:
                    sids[tid] = []
                    new_subsets[tid] = []
                sids[tid].append(subset.id)
                new_subsets[tid].append(subset)

        if create_subsets:
            return new_templates, new_subsets
        else:
            return new_templates

    def get_molecule_smiles(self):
        """Return the a list of the canonical SMILES for each molecule..

        Returns
        -------
        [str]
            The canonical SMILES for each molecule, in order that they are found.
        """
        # Find the molecules
        molecules = self.find_molecules()
        n_molecules = len(molecules)

        # And the molecule each atom is in
        atom_to_molecule = {}
        for molecule, atoms in enumerate(molecules):
            for atom in atoms:
                atom_to_molecule[atom] = molecule

        # The bonds in each molecule
        bonds_per_molecule = [[] for i in range(n_molecules)]
        for bond in self.bonds.bonds():
            i = bond["i"]
            j = bond["j"]
            order = bond["bondorder"]
            molecule = atom_to_molecule[i]
            bonds_per_molecule[molecule].append((i, j, order))

        # Get the canonical smiles for each molecule
        to_can = openbabel.OBConversion()
        to_can.SetOutFormat("can")
        ob_mol = openbabel.OBMol()
        atnos = self.atoms.atomic_numbers

        atom_index = {j: i for i, j in enumerate(self.atoms.ids)}

        SMILES = []
        for molecule, atoms in enumerate(molecules):
            to_index = {j: i for i, j in enumerate(atoms)}
            molecule_atnos = [atnos[atom_index[i]] for i in atoms]

            ob_mol.Clear()
            for atom, atno in zip(atoms, molecule_atnos):
                ob_atom = ob_mol.NewAtom()
                ob_atom.SetAtomicNum(atno)
            bonds = []
            for i, j, order in bonds_per_molecule[molecule]:
                bonds.append((to_index[i], to_index[j], order))
                # 1-based indices in ob.
                ob_mol.AddBond(to_index[i] + 1, to_index[j] + 1, order)

            canonical = to_can.WriteString(ob_mol).strip()
            SMILES.append(canonical)
        return SMILES

    def reimage_bonded_atoms(self, reimage_molecules=True):
        """Ensure that the atoms in a molecule are "near" each other.

        In a periodic system atoms can be translated by a unit cell without changing
        anything. However the bonds in a molecule need to account for such translations
        of atoms otherwise a bond may appear to be very long.

        It is often convenient physically move atoms by the correct cell translations
        to bring the atoms of a molecule close to each other. That is what this method
        does, moving all the atoms close to the first. Optionally the molecule is
        also moved to bring its geometric center into the primary unit cell, i.e. with
        fractional coordinates in the range [0..1).

        Parameters
        ----------
        reimage_molecules : bool = True
            Whether to move molecules into the primary unit cell.

        Returns
        -------
        bool
            True if the coordinates were changed.
        """
        if self.periodicity == 0:
            # Nothing to do.
            return False

        logger.log(0, f"Configuration {self.id} {self}")
        logger.log(0, self.atoms)
        logger.log(0, 10 * "+")
        # The bonds from each atom
        bonds = self.bonded_neighbors(as_indices=True)

        logger.debug("Bonds:")
        logger.debug(pprint.pformat(bonds))
        logger.debug("")

        # And coordinates as fractionals
        tmp = self.atoms.get_coordinates(fractionals=True)
        xyz = [[x, y, z] for x, y, z in tmp]

        visited = [False] * len(xyz)

        # Find the molecules
        molecules = self.find_molecules(as_indices=True)

        # Work through the atoms, see if bonded atoms need to be moved
        changed = False
        for atoms in molecules:
            logger.debug("")
            logger.debug(f"{atoms=}")
            # if len(atoms) > 1:
            #     atom = sorted(atoms)[0]
            #     logger.debug(f" --> {atom}")
            #     if self._image_helper(bonds, xyz, atom):
            #         changed = True

            changed = self._image_helper2(bonds, xyz, visited, atoms[0], changed)
        # Check temporarily.
        for iatom, jatoms in enumerate(bonds):
            xyzi = xyz[iatom]
            for jatom in jatoms:
                xyzj = xyz[jatom]
                for vi, vj in zip(xyzi, xyzj):
                    delta = floor(vj - vi + 0.5)
                    if delta != 0:
                        print(
                            f"Warning: {iatom} - {jatom} delta = {vj:.3f} - {vi:.3f} = "
                            f"{delta}"
                        )

        # Move the center of molecules into the primary cell, if requested.
        if reimage_molecules:
            for atoms in molecules:
                cx = cy = cz = 0.0
                for atom in atoms:
                    x, y, z = xyz[atom]
                    cx += x
                    cy += y
                    cz += z
                n = len(atoms)
                dx = int(cx / n)
                dy = int(cy / n)
                dz = int(cz / n)
                if dx != 0 or dy != 0 or dz != 0:
                    changed = True
                    for atom in atoms:
                        x, y, z = xyz[atom]
                        xyz[atom] = [x - dx, y - dy, z - dz]
        if changed:
            self.atoms.set_coordinates(xyz, fractionals=True)

        return changed

    def _image_helper(self, bonds, xyz, atom):
        """Translate all the bonded neighbors with larger index.

        Parameters
        ----------
        bonds : [[int]]
            The bonded neighbors for each atom.
        xyz : [[float * 3]]
            The fractional coordinates of the atoms.
        """
        logger.debug(f"       helper {atom}: {bonds[atom]}")
        changed = False
        xyzi = xyz[atom]
        for jatom in bonds[atom]:
            logger.debug(f"\t\t{atom=} {jatom}")
            if jatom < atom:
                continue
            xyzj = xyz[jatom]
            new_xyz = []
            shift = False
            for vi, vj in zip(xyzi, xyzj):
                delta = floor(vj - vi + 0.5)
                if delta == 0:
                    new_xyz.append(vj)
                else:
                    new_xyz.append(vj - delta)
                    shift = True
                logger.log(
                    0, f"{atom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
            if shift:
                changed = True
                logger.debug(
                    f"{atom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
                xyz[jatom] = new_xyz

        for jatom in bonds[atom]:
            if jatom < atom:
                continue
            if self._image_helper(bonds, xyz, jatom):
                changed = True

        return changed

    def reimage_molecules(self):
        """Reimage molecules into the primary unit cell.

        The molecules are moved to bring their geometric center into the primary unit
        cell, i.e. with fractional coordinates in the range [0..1).

        Returns
        -------
        bool
            True if the coordinates were changed.
        """
        if self.periodicity == 0:
            # Nothing to do.
            return False

        # And coordinates as fractionals
        xyz = self.atoms.get_coordinates(fractionals=True)

        # Find the molecules
        molecules = self.find_molecules(as_indices=True)

        # Move the centor of molecules into the primary cell, if requested.
        changed = False
        for atoms in molecules:
            cx = cy = cz = 0.0
            for atom in atoms:
                x, y, z = xyz[atom]
                cx += x
                cy += y
                cz += z
            n = len(atoms)
            dx = int(cx / n)
            dy = int(cy / n)
            dz = int(cz / n)
            if dx != 0 or dy != 0 or dz != 0:
                changed = True
                for atom in atoms:
                    x, y, z = xyz[atom]
                    xyz[atom] = [x - dx, y - dy, z - dz]
        if changed:
            self.atoms.set_coordinates(xyz, fractionals=True)

        return changed

    def _image_helper2(self, bonds, xyz, visited, iatom, changed):
        visited[iatom] = True
        xyzi = xyz[iatom]
        for jatom in bonds[iatom]:
            logger.debug(f"\t\t{iatom=} {jatom}")
            if visited[jatom]:
                continue
            new_xyz = []
            xyzj = xyz[jatom]
            shift = False
            for vi, vj in zip(xyzi, xyzj):
                delta = floor(vj - vi + 0.5)
                if delta == 0:
                    new_xyz.append(vj)
                else:
                    new_xyz.append(vj - delta)
                    changed = True
                    shift = True
            if shift:
                logger.debug(
                    f"{iatom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
            xyz[jatom] = new_xyz
            changed = self._image_helper2(bonds, xyz, visited, jatom, changed)
        return changed

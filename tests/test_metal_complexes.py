# -*- coding: utf-8 -*-

"""Metal complexes, ionic + dative (Paul, 2026-11-19): every metal-ligand bond is a
dative bond from the ligand's donor atom to the metal (MDL bond type 9), each
ligand has its closed-shell Lewis structure with the least charge separation, the
metals carry the rest of the net charge, and metal-metal bonds are single bonds.

Dative bonds keep their direction, although bonds are stored with i < j (bond
orders 9 and 10); before, Fe <- C became Fe -> C on writing, and RDKit refused the
file. The geometries here are idealized, built for the tests."""

import math

import numpy as np
import pytest

from molsystem import DATIVE, DATIVE_REVERSED, dative_donor_acceptor

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402


def _build(configuration, symbols, xyz, charge=0, multiplicity=1):
    xyz = np.asarray(xyz, dtype=float)
    configuration.atoms.append(
        x=xyz[:, 0].tolist(), y=xyz[:, 1].tolist(), z=xyz[:, 2].tolist(), symbol=symbols
    )
    configuration.charge = charge
    configuration.spin_multiplicity = multiplicity
    configuration.perceive_bonds()
    configuration.assign_bond_orders(charge=charge)
    return configuration


def _smiles(configuration):
    return Chem.MolToSmiles(Chem.RemoveHs(configuration.to_RDKMol()))


def _iron_pentacarbonyl():
    """Trigonal bipyramidal Fe(CO)5: Fe-C 1.81, C-O 1.15 Å."""
    symbols = ["Fe"]
    xyz = [[0.0, 0.0, 0.0]]
    directions = [(0, 0, 1), (0, 0, -1)] + [
        (math.cos(a), math.sin(a), 0) for a in (0, 2 * math.pi / 3, 4 * math.pi / 3)
    ]
    for d in directions:
        d = np.array(d, dtype=float)
        symbols += ["C", "O"]
        xyz += [1.81 * d, (1.81 + 1.15) * d]
    return symbols, xyz


def _ferrocene():
    """Eclipsed ferrocene: Cp rings 1.66 Å from Fe, C-C 1.42, C-H 1.08 Å."""
    r = 1.42 / (2 * math.sin(math.pi / 5))
    symbols = ["Fe"]
    xyz = [[0.0, 0.0, 0.0]]
    for z in (1.66, -1.66):
        for k in range(5):
            a = 2 * math.pi * k / 5
            symbols.append("C")
            xyz.append([r * math.cos(a), r * math.sin(a), z])
        for k in range(5):
            a = 2 * math.pi * k / 5
            symbols.append("H")
            xyz.append([(r + 1.08) * math.cos(a), (r + 1.08) * math.sin(a), z])
    return symbols, xyz


def test_the_direction_of_a_dative_bond_is_kept(configuration):
    """A dative bond given from the higher-numbered atom is stored as 10."""
    configuration.atoms.append(
        x=[0.0, 1.81, 2.96],
        y=[0.0, 0.0, 0.0],
        z=[0.0, 0.0, 0.0],
        symbol=["Fe", "C", "O"],
    )
    fe, c, o = configuration.atoms.ids
    configuration.bonds.append(i=[c, c], j=[fe, o], bondorder=[DATIVE, 3])
    rows = {(r["i"], r["j"]): r["bondorder"] for r in configuration.bonds.bonds()}
    assert rows == {(fe, c): DATIVE_REVERSED, (c, o): 3}
    assert dative_donor_acceptor(fe, c, DATIVE_REVERSED) == (c, fe)
    assert dative_donor_acceptor(c, o, 3) is None


def test_dative_bonds_survive_sdf_and_rdkit(configuration, system):
    symbols, xyz = _iron_pentacarbonyl()
    _build(configuration, symbols, xyz)
    expected = "[O+]#[C-]->[Fe](<-[C-]#[O+])(<-[C-]#[O+])(<-[C-]#[O+])<-[C-]#[O+]"
    assert _smiles(configuration) == expected
    # MDL type 9 from C (donor) to Fe: RDKit reads the written file
    text = configuration.to_sdf_text()
    mol = Chem.MolFromMolBlock(text, removeHs=False)
    assert mol is not None
    assert Chem.MolToSmiles(mol, isomericSmiles=False) == expected
    # and SEAMM reads it back the same way
    other = system.create_configuration(name="again")
    other.from_sdf_text(text)
    assert _smiles(other) == expected
    # and from RDKit
    third = system.create_configuration(name="from RDKit")
    third.from_RDKMol(configuration.to_RDKMol())
    assert _smiles(third) == expected


def test_iron_pentacarbonyl(configuration):
    symbols, xyz = _iron_pentacarbonyl()
    _build(configuration, symbols, xyz)
    orders = sorted(r["bondorder"] for r in configuration.bonds.bonds())
    assert orders.count(3) == 5 and sum(o in (9, 10) for o in orders) == 5
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges[0] == 0  # Fe(0)
    assert sorted(charges[1:]) == [-1] * 5 + [1] * 5  # [C-]#[O+]


def test_ferrocene(configuration):
    """Two Cp- (aromatic), each C bonded to Fe by a dative bond, and Fe2+."""
    symbols, xyz = _ferrocene()
    _build(configuration, symbols, xyz)
    datives = [r for r in configuration.bonds.bonds() if r["bondorder"] in (9, 10)]
    assert len(datives) == 10
    for row in datives:
        donor, acceptor = dative_donor_acceptor(row["i"], row["j"], row["bondorder"])
        assert configuration.atoms.symbols[acceptor - 1] == "Fe"
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges[0] == 2 and sum(charges) == 0
    mol = configuration.to_RDKMol()
    assert Chem.MolToSmiles(mol).count("->") + Chem.MolToSmiles(mol).count("<-") == 10


def test_open_shell_complex_and_halides(configuration, system):
    """Quintet Cr(OH)2: OH- ligands, Cr2+; the unpaired electrons are the metal's,
    so no radicals on the ligands. TiCl4: Cl- and Ti4+."""
    _build(
        configuration,
        ["Cr", "O", "H", "O", "H"],
        [[0, 0, 0], [1.85, 0, 0], [2.30, 0.88, 0], [-1.85, 0, 0], [-2.30, -0.88, 0]],
        multiplicity=5,
    )
    assert _smiles(configuration) == "[OH-]->[Cr+2]<-[OH-]"
    mol = configuration.to_RDKMol()
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 0
    assert mol.GetIntProp("SEAMM|spin multiplicity|int|") == 5

    titanium = system.create_configuration(name="TiCl4")
    d = 2.17 / math.sqrt(3)
    _build(
        titanium,
        ["Ti", "Cl", "Cl", "Cl", "Cl"],
        [[0, 0, 0], [d, d, d], [d, -d, -d], [-d, d, -d], [-d, -d, d]],
    )
    assert _smiles(titanium) == "[Cl-]->[Ti+4](<-[Cl-])(<-[Cl-])<-[Cl-]"


def test_metal_metal_bonds_are_single(configuration):
    """Hg2Cl2: Cl- -> Hg+ - Hg+ <- Cl-."""
    _build(
        configuration,
        ["Cl", "Hg", "Hg", "Cl"],
        [[-3.515, 0, 0], [-1.265, 0, 0], [1.265, 0, 0], [3.515, 0, 0]],
    )
    rows = sorted(
        (
            configuration.atoms.symbols[r["i"] - 1],
            configuration.atoms.symbols[r["j"] - 1],
            r["bondorder"],
        )
        for r in configuration.bonds.bonds()
    )
    assert ("Hg", "Hg", 1) in rows
    assert _smiles(configuration) == "[Cl-]->[Hg+][Hg+]<-[Cl-]"


def test_charged_complex(configuration):
    """[FeCl4]-: four Cl- and Fe3+ (the metal takes the rest of the charge)."""
    d = 2.19 / math.sqrt(3)
    _build(
        configuration,
        ["Fe", "Cl", "Cl", "Cl", "Cl"],
        [[0, 0, 0], [d, d, d], [d, -d, -d], [-d, d, -d], [-d, -d, d]],
        charge=-1,
        multiplicity=6,
    )
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges == [3, -1, -1, -1, -1]


def test_the_center_of_a_chelate_ring_is_no_donor(configuration):
    """[La(kappa2-HCOO)]2+: the wider metal tolerance reaches the formate C, 0.34 Å
    further than the two O; that bond is dropped, so O,O' are the donors."""
    _build(
        configuration,
        ["La", "O", "O", "C", "H"],
        [[0, 0, 0], [2.30, 1.10, 0], [2.30, -1.10, 0], [2.894, 0, 0], [3.984, 0, 0]],
        charge=2,
    )
    pairs = {
        tuple(
            sorted(
                (
                    configuration.atoms.symbols[r["i"] - 1],
                    configuration.atoms.symbols[r["j"] - 1],
                )
            )
        )
        for r in configuration.bonds.bonds()
    }
    assert ("C", "La") not in pairs and ("La", "O") in pairs
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges[0] == 3 and sum(charges) == 2


def test_lone_atom_ligands_are_anions(configuration, system):
    """A lone H, C, N on a metal is H-, C4-, N3- (PM-8: FeH, ZrC, ZrN), not a
    covalent hydride or a bare atom."""
    _build(configuration, ["Fe", "H"], [[0, 0, 0], [1.6, 0, 0]], multiplicity=4)
    assert _smiles(configuration) == "[H-]->[Fe+]"
    carbide = system.create_configuration(name="ZrC")
    _build(carbide, ["Zr", "C"], [[0, 0, 0], [1.8, 0, 0]], multiplicity=3)
    assert carbide.atoms.get_column_data("formal_charge") == [4, -4]
    nitride = system.create_configuration(name="ZrN")
    _build(nitride, ["Zr", "N"], [[0, 0, 0], [1.7, 0, 0]], multiplicity=2)
    assert nitride.atoms.get_column_data("formal_charge") == [3, -3]


def test_alkali_metals_take_their_charge_first(configuration):
    """Na bonded in NaZrF5 (when the alkali metals are not left out of the
    perception) is Na+, and Zr takes the rest: Zr4+, not Na2+/Zr3+."""
    d = 1.95
    symbols = ["Zr", "F", "F", "F", "F", "F", "Na"]
    xyz = [[0, 0, 0], [d, 0, 0], [-d, 0, 0], [0, d, 0], [0, -d, 0], [0, 0, d]]
    xyz.append([0, 0, d + 2.2])
    xyz = np.asarray(xyz, dtype=float)
    configuration.atoms.append(
        x=xyz[:, 0].tolist(), y=xyz[:, 1].tolist(), z=xyz[:, 2].tolist(), symbol=symbols
    )
    configuration.charge = 0
    configuration.spin_multiplicity = 1
    configuration.perceive_bonds(exclude=())
    configuration.assign_bond_orders(charge=0)
    charges = configuration.atoms.get_column_data("formal_charge")
    assert charges[0] == 4 and charges[-1] == 1 and sum(charges) == 0

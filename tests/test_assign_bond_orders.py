# -*- coding: utf-8 -*-

"""Tests of TopologyMixin.assign_bond_orders: perceived bonds are single bonds, so a
structure read from a file without bond orders (extended XYZ) needs its bond orders
assigned, or ethylene reads as [CH2+][CH2-] and the Thermochemistry step builds the
wrong CBH fragments (thermo, ChemAI 5585, 2026-10-31)."""

import pytest


def _connectivity_only(configuration, smiles, charge=0, multiplicity=1):
    """The molecule's geometry with perceived single bonds, as from extended XYZ."""
    configuration.from_smiles(smiles)
    configuration.bonds.delete()
    if "formal_charge" in configuration.atoms:
        configuration.atoms["formal_charge"][0:] = [0] * configuration.atoms.n_atoms
    configuration.charge = charge
    configuration.spin_multiplicity = multiplicity
    configuration.perceive_bonds()
    return configuration


@pytest.mark.parametrize(
    "smiles, charge, before, after",
    [
        ("C=C", 0, "[CH2+][CH2-]", "C=C"),
        ("CC(C)=O", 0, "C[C+](C)[O-]", "CC(C)=O"),
        ("Cc1ccccc1", 0, None, "Cc1ccccc1"),
        ("CC(=O)[O-]", -1, None, "CC(=O)[O-]"),
        ("O=C1OCCO1.O=C1OCCO1", 0, None, "O=C1OCCO1.O=C1OCCO1"),
    ],
)
def test_bond_orders_are_assigned(configuration, smiles, charge, before, after):
    _connectivity_only(configuration, smiles, charge)
    if before is not None:
        assert configuration.canonical_smiles == before
    n_ok, n_failed = configuration.assign_bond_orders(charge=charge)
    assert n_failed == 0 and n_ok == smiles.count(".") + 1
    assert configuration.canonical_smiles == after


def test_ions_and_hypervalent_molecules(configuration):
    """Li+ is an ion; PF6- (hypervalent) cannot be assigned by RDKit and keeps its
    single bonds, which are its Lewis structure."""
    _connectivity_only(configuration, "[Li+].F[P-](F)(F)(F)(F)F", charge=0)
    n_ok, n_failed = configuration.assign_bond_orders(charge=0)
    assert (n_ok, n_failed) == (1, 1)
    assert configuration.canonical_smiles == "F[P-](F)(F)(F)(F)F.[Li+]"


def test_open_shell_molecules_and_periodic_configurations(configuration):
    """Radicals get their Lewis structures too (the methyl radical, with one
    radical electron); periodic configurations are left alone."""
    _connectivity_only(configuration, "[CH3]", multiplicity=2)
    assert configuration.assign_bond_orders() == (1, 0)
    assert configuration.atoms.get_column_data("radical_electrons")[0] == 1
    assert configuration.canonical_smiles == "[CH3]"
    configuration.periodicity = 3
    configuration.spin_multiplicity = 1
    assert configuration.assign_bond_orders() == (0, 0)


def test_the_net_charge_of_a_single_molecule_is_kept(configuration):
    """Acetate declared neutral is the acetoxyl radical, not the anion."""
    _connectivity_only(configuration, "CC(=O)[O-]", charge=0)
    assert configuration.assign_bond_orders(charge=0) == (1, 0)
    assert sum(configuration.atoms.get_column_data("formal_charge")) == 0


def test_charges_that_do_not_add_up_are_reported(configuration, caplog):
    _connectivity_only(configuration, "CC(=O)[O-].C", charge=0)
    configuration.assign_bond_orders(charge=0)
    assert "add up to -1, not the net charge 0" in caplog.text


def _atoms(configuration, symbols, xyz, charge=0, multiplicity=1):
    configuration.atoms.append(
        x=[p[0] for p in xyz],
        y=[p[1] for p in xyz],
        z=[p[2] for p in xyz],
        symbol=symbols,
    )
    configuration.charge = charge
    configuration.spin_multiplicity = multiplicity
    configuration.perceive_bonds()
    configuration.assign_bond_orders(charge=charge)
    return configuration


def _rdkit_smiles(configuration):
    from rdkit import Chem

    return Chem.MolToSmiles(configuration.to_RDKMol())


@pytest.mark.parametrize(
    "symbols, xyz, multiplicity, smiles",
    [
        # no closed-shell Lewis structure: single bonds and radicals, not
        # [Ge+][Ge-] or [C+2][Ge-2]
        (["Ge", "Ge"], [[0, 0, 0], [2.44, 0, 0]], 1, "[Ge][Ge]"),
        (["Ge", "C"], [[0, 0, 0], [1.85, 0, 0]], 1, "[C][Ge]"),
        (
            ["Ge", "Ge", "Ge"],
            [[0, 0, 0], [2.4, 0, 0], [1.2, 2.08, 0]],
            1,
            "[Ge]1[Ge][Ge]1",
        ),
        # a closed-shell structure is kept whatever the multiplicity (triplet
        # AsN, SeTe), rather than charges fitted to it ([N+]#[As-])
        (["As", "N"], [[0, 0, 0], [1.62, 0, 0]], 3, "N#[As]"),
        (["Se", "Te"], [[0, 0, 0], [2.4, 0, 0]], 3, "[Se]=[Te]"),
        (
            ["Te", "Ge", "Te"],
            [[-2.3, 0, 0], [0, 0, 0], [2.3, 0, 0]],
            3,
            "[Te]=[Ge]=[Te]",
        ),
    ],
)
def test_species_without_a_simple_lewis_structure(
    configuration, symbols, xyz, multiplicity, smiles
):
    _atoms(configuration, symbols, xyz, multiplicity=multiplicity)
    assert _rdkit_smiles(configuration) == smiles


def test_radicals_survive_sdf_and_rdkit(configuration, system):
    _atoms(configuration, ["Ge", "Ge"], [[0, 0, 0], [2.44, 0, 0]])
    text = configuration.to_sdf_text()
    other = system.create_configuration(name="from SDF")
    other.from_sdf_text(text)
    assert other.atoms.get_column_data("radical_electrons") == [3, 3]
    assert _rdkit_smiles(other) == "[Ge][Ge]"
    third = system.create_configuration(name="from RDKit")
    third.from_RDKMol(configuration.to_RDKMol())
    assert _rdkit_smiles(third) == "[Ge][Ge]"


def test_vinyl_and_ethyl_radicals(configuration, system):
    """Radicals with a structure that has the right number of radical electrons
    keep it: C=[CH], not a triradical."""
    _connectivity_only(configuration, "C=[CH]", multiplicity=2)
    configuration.assign_bond_orders(charge=0)
    assert _rdkit_smiles(configuration) == "[H][C]=C([H])[H]"
    ethyl = system.create_configuration(name="ethyl")
    _connectivity_only(ethyl, "[CH2]C", multiplicity=2)
    ethyl.assign_bond_orders(charge=0)
    assert _rdkit_smiles(ethyl) == "[H][C]([H])C([H])([H])[H]"

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


def test_open_shell_and_periodic_are_left_alone(configuration):
    _connectivity_only(configuration, "[CH3]", multiplicity=2)
    assert configuration.assign_bond_orders() == (0, 0)
    configuration.periodicity = 3
    configuration.spin_multiplicity = 1
    assert configuration.assign_bond_orders() == (0, 0)


def test_charges_that_do_not_add_up_are_reported(configuration, caplog):
    _connectivity_only(configuration, "CC(=O)[O-]", charge=0)
    configuration.assign_bond_orders(charge=0)
    assert "add up to -1, not the net charge 0" in caplog.text

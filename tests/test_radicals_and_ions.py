#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""SMILES of radicals and ions: the charges and unpaired electrons of the
configuration carried into RDKit and back out, from the atoms' stored formal
charges or, without them, from the net charge and spin multiplicity."""

import pytest

from rdkit import Chem

CASES = [
    # (SMILES, charge, spin multiplicity)
    ("[CH2]C", 0, 2),
    ("C[O]", 0, 2),
    ("[CH3]", 0, 2),
    ("[CH2]", 0, 3),
    ("[OH]", 0, 2),
    ("[NH2]", 0, 2),
    ("C[CH]C", 0, 2),
    ("[O]O", 0, 2),
    ("[O][O]", 0, 3),
    ("[CH2][CH2]", 0, 3),
    ("[N]=O", 0, 2),
    ("[CH2]C=C", 0, 2),
    ("C=CC[O]", 0, 2),
    ("[NH3+]C", 1, 1),
    ("CC(=O)[O-]", -1, 1),
    ("[CH3+]", 1, 1),
    ("[CH2+]C", 1, 1),
    ("C[NH2+]", 1, 2),
    ("[CH3-]", -1, 1),
    ("C[O-]", -1, 1),
    ("[NH4+]", 1, 1),
    ("[BH4-]", -1, 1),
    ("C[N+](=O)[O-]", 0, 1),
    ("C[S+](C)[O-]", 0, 1),
    ("[C-]#[O+]", 0, 1),
    ("[NH3+]CC(=O)[O-]", 0, 1),
    ("OS(=O)(=O)O", 0, 1),
    ("CS(C)=O", 0, 1),
    ("F[P-](F)(F)(F)(F)F", -1, 1),
    ("c1ccccc1", 0, 1),
    ("[c]1ccccc1", 0, 2),
    ("[CH2]c1ccccc1", 0, 2),
    ("c1cc[nH]c1", 0, 1),
    ("c1cc[nH+]cc1", 1, 1),
    ("[cH-]1cccc1", -1, 1),
    ("[O-][n+]1ccccc1", 0, 1),
    ("O=c1cc[nH]cc1", 0, 1),
    ("Cc1ccccc1[O-]", -1, 1),
]


def canonical(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def without_formal_charges(db, configuration):
    """A copy of the configuration with only the atoms, coordinates and bonds, as
    from a file that has no formal charges."""
    other = db.create_system().create_configuration()
    xyz = configuration.atoms.get_coordinates()
    ids = other.atoms.append(
        x=[p[0] for p in xyz],
        y=[p[1] for p in xyz],
        z=[p[2] for p in xyz],
        atno=configuration.atoms.atomic_numbers,
    )
    index = dict(zip(configuration.atoms.ids, ids))
    bonds = list(configuration.bonds.bonds())
    other.bonds.append(
        i=[index[row["i"]] for row in bonds],
        j=[index[row["j"]] for row in bonds],
        bondorder=[row["bondorder"] for row in bonds],
    )
    other.charge = configuration.charge
    other.spin_multiplicity = configuration.spin_multiplicity
    return other


@pytest.mark.parametrize("smiles, charge, multiplicity", CASES)
def test_round_trip(configuration, smiles, charge, multiplicity):
    """from_smiles sets the charge and multiplicity; to_smiles gives the SMILES."""
    configuration.from_smiles(smiles)
    assert configuration.charge == charge
    assert configuration.spin_multiplicity == multiplicity
    assert configuration.to_smiles(canonical=True, isomeric=True) == canonical(smiles)


@pytest.mark.parametrize("smiles, charge, multiplicity", CASES)
def test_from_charge_and_spin(db, configuration, smiles, charge, multiplicity):
    """Without formal charges on the atoms they follow from the net charge and spin."""
    configuration.from_smiles(smiles)
    other = without_formal_charges(db, configuration)
    assert other.to_smiles(canonical=True, isomeric=True) == canonical(smiles)


def test_charges_follow_the_configuration(configuration):
    """The same atoms and bonds are a radical, a cation or an anion."""
    configuration.from_smiles("[CH2]C")
    for charge, multiplicity, smiles in (
        (0, 2, "[CH2]C"),
        (1, 1, "[CH2+]C"),
        (-1, 1, "[CH2-]C"),
    ):
        configuration.charge = charge
        configuration.spin_multiplicity = multiplicity
        assert configuration.to_smiles(canonical=True) == smiles


def test_rdkit_molecule(configuration):
    """The RDKit molecule has the radical electrons and formal charges."""
    configuration.from_smiles("C[N+](=O)[O-]")
    mol = configuration.to_RDKMol()
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 0
    charges = {
        a.GetSymbol(): sorted(
            b.GetFormalCharge()
            for b in mol.GetAtoms()
            if b.GetSymbol() == a.GetSymbol()
        )
        for a in mol.GetAtoms()
    }
    assert charges == {"C": [0], "N": [1], "O": [-1, 0], "H": [0, 0, 0]}
    assert Chem.GetFormalCharge(mol) == 0

    configuration.from_smiles("[CH2]C")
    mol = configuration.to_RDKMol()
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 1
    assert Chem.GetFormalCharge(mol) == 0


def test_electrolyte(db, configuration):
    """Several ions and a radical in one structure, charges inferred."""
    smiles = ".".join(
        ["O"] * 20
        + ["[Li+]"] * 3
        + ["F[P-](F)(F)(F)(F)F", "CC(=O)[O-]", "[NH3+]CC(=O)[O-]", "[CH2]C"]
    )
    configuration.from_smiles(smiles)
    assert configuration.charge == 1
    assert configuration.spin_multiplicity == 2
    expected = configuration.to_smiles(canonical=True)
    assert expected == canonical(smiles)
    other = without_formal_charges(db, configuration)
    assert other.to_smiles(canonical=True) == expected

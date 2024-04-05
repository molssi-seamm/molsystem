#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

import molsystem  # noqa: F401

"""Tests for templates in the `molsystem` package."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]


def test_construction(db):
    """Simplest test that we can make a new templates object"""
    templates = db.templates
    assert str(type(templates)) == "<class 'molsystem.templates._Templates'>"


def test_create(db):
    """Test that we can make a new template"""
    templates = db.templates
    template = templates.create(name="H2O", category="molecule")
    assert template.id == 1


def test_create_many(db):
    """Test that we can make several templates."""
    templates = db.templates
    names = ["H2O", "H2", "O2"]
    new = templates.create_many(name=names, category="molecule")
    assert [x.id for x in new] == [1, 2, 3]
    assert templates.categories == ["molecule"]
    assert templates.names("molecule") == sorted(names)


def test_making_full_templates(amino_acids):
    """Create full templates."""
    answer = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
    templates = amino_acids.templates
    for system in amino_acids.systems:
        templates.create(
            name=system.name,
            category="amino acid",
            configuration=system.configuration.id,
        )
    assert amino_acids.n_templates == 20
    assert templates.names("amino acid") == answer


def test_get_template(aa_templates):
    """Test getting a template by index and name."""
    template = aa_templates.get(8)
    assert template.name == "GLY"
    template = aa_templates.get("GLU", category="amino acid")
    assert template.id == 6


def test_template(gly):
    """Test the gly template."""
    assert gly.n_atoms == 10
    assert gly.n_bonds == 9


def test_template_coordinate_system(gly):
    """Test the gly template coordinate system."""
    assert gly.coordinate_system == "Cartesian"


def test_template_formula(gly):
    """Test the gly template formula."""
    assert gly.formula == ("C2 H5 N O2", "C2 H5 N O2", 1)


def test_template_mass(gly):
    """Test the gly template mass."""
    assert abs(gly.mass - 75.067) <= 0.001


def test_template_periodicity(gly):
    """Test the gly template periodicity."""
    assert gly.periodicity == 0


def test_template_atoms(gly):
    """Test the atoms for the template."""
    answer = ["N", "CA", "C", "O", "OXT", "H", "H2", "HA2", "HA3", "HXT"]
    names = gly.atoms.get_column_data("name")
    if names != answer:
        print(names)
    assert names == answer


def test_atoms_data(gly):
    """Test the data for the atoms for the template."""
    answer = {
        "id": [127, 128, 129, 130, 131, 132, 133, 134, 135, 136],
        "atno": [7, 6, 6, 8, 8, 1, 1, 1, 1, 1],
        "name": ["N", "CA", "C", "O", "OXT", "H", "H2", "HA2", "HA3", "HXT"],
        "formal_charge": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "_chem_comp_atom.alt_atom_id": [
            "N",
            "CA",
            "C",
            "O",
            "OXT",
            "H",
            "HN2",
            "HA1",
            "HA2",
            "HXT",
        ],
        "_chem_comp_atom.pdbx_align": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "_chem_comp_atom.pdbx_aromatic_flag": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "_chem_comp_atom.pdbx_leaving_atom_flag": [1, 1, 1, 1, 0, 1, 0, 1, 1, 0],
        "_chem_comp_atom.pdbx_stereo_config": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "_chem_comp_atom.pdbx_component_atom_id": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "configuration": [8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
        "x": [
            25.463,
            25.329,
            26.081,
            27.024,
            25.702,
            25.494,
            26.307,
            24.27,
            25.731,
            26.236,
        ],
        "y": [
            35.609,
            37.024,
            37.335,
            36.627,
            38.256,
            35.15,
            35.421,
            37.305,
            37.59,
            38.3,
        ],
        "z": [
            47.047,
            46.85,
            45.572,
            45.222,
            44.874,
            46.159,
            47.549,
            46.757,
            47.703,
            44.09,
        ],
    }
    for column in gly.atoms:
        if column in ("gx", "gy", "gz", "vx", "vy", "vz"):
            continue
        result = gly.atoms.get_column_data(column)
        if column not in answer:
            print(f"'{column}': {result}")
        assert column in answer
        if result != answer[column]:
            print(f"'{column}': {result}")
        assert result == answer[column]


def test_template_bonds(gly):
    """Test the bonds for the template."""
    answer = [1, 1, 1, 1, 1, 1, 2, 1, 1]
    bondorders = gly.bonds.get_column_data("bondorder")
    if bondorders != answer:
        print(bondorders)
    assert bondorders == answer


def test_print_template_bonds(gly):
    """Test the string reprentation of the bonds in the template."""
    answer = """\
       i    j  bondorder symop1 symop2 _chem_comp_bond.comp_id
120  127  128          1      .      .                     GLY
121  127  132          1      .      .                     GLY
122  127  133          1      .      .                     GLY
123  128  129          1      .      .                     GLY
124  128  134          1      .      .                     GLY
125  128  135          1      .      .                     GLY
126  129  130          2      .      .                     GLY
127  129  131          1      .      .                     GLY
128  131  136          1      .      .                     GLY"""

    result = gly.bonds.to_dataframe().to_string()
    if result != answer:
        print(result)
    assert result == answer

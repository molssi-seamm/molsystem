#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

import molsystem  # noqa: F401
"""Tests for `molsystem` package."""

x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [7.0, 8.0, 9.0]
atno = [8, 1, 1]


def test_construction(system):
    """Simplest test that we can make a new template"""
    templates = system['template']
    templates.append(name='H2O', type='molecule')
    assert str(type(templates)) == "<class 'molsystem.template._Template'>"


def test_set_current_template(system):
    """Test that we can set the current template."""
    templates = system['template']
    tid = templates.append(name='H2O', type='molecule')[0]
    templates.current_template = tid
    assert templates.n_rows == 2 and templates.current_template == tid


def test_adding_atoms(system):
    """Test that we can add atoms to a new template."""
    templates = system['template']
    tid = templates.append(name='H2O', type='molecule')[0]
    templates.current_template = tid

    atoms = system['templateatom']
    atoms.append(atno=atno, x=x, y=y, z=z)
    assert atoms.n_atoms == 3


def test_template_str(templates):
    """Print the water template."""
    atom_string = """\
   template name  atno        x    y         z
1         2    O     8  0.00000  0.0  0.000000
2         2   H1     1  0.75695  0.0  0.585882
3         2   H2     1 -0.75695  0.0  0.585882"""

    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    atoms = system['templateatoms']
    if str(atoms) != atom_string:
        print(atoms)
    assert str(atoms) == atom_string


def test_template_n_atoms(templates):
    """Get the number of atoms in acetic acid template."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    atoms = system['templateatoms']
    assert atoms.n_atoms == 8


def test_invalid_key(templates):
    """Test the error trying to get a non-existant template."""
    system = templates
    templates = system['template']
    with pytest.raises(KeyError) as e:
        templates.current_template = 99
    assert str(e.value) == '''"Template '99' does not exist."'''


def test_invalid_key_2(templates):
    """Test the error trying to get a non-existant template."""
    answer = '''"There is no template 'junk' of type 'molecule'."'''
    system = templates
    templates = system['template']
    with pytest.raises(KeyError) as e:
        templates.set_current_template('junk', type_='molecule')
    assert str(e.value) == answer


def test_get_column_x(templates):
    """Test getting a column from the atoms in a template."""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    atoms = system['templateatom']
    x = atoms['x']
    assert x.equal([0.0, 0.75695, -0.75695], 1.0e-05)


def test_get_column_atno(templates):
    """Test getting a column from the atoms in a template."""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    atoms = system['templateatom']
    atno = atoms['atno']
    assert atno == [8, 1, 1]


def test_get_column_error(templates):
    """Test getting a column that doesn't exist."""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    atoms = system['templateatom']
    with pytest.raises(KeyError) as e:
        junk = atoms['junk']  # noqa: F841
    assert str(e.value) == '''"'junk' not in template atoms"'''


def test_select_atoms(templates):
    """Select just some atoms in the template."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    atoms = system['templateatom']
    names = []
    x = []
    for row in atoms.atoms('atno', '==', 8):
        names.append(row['name'])
        x.append(row['x'])
    assert names == ['O', 'OH']
    assert x == [-0.1323, 0.9757]


def test_select_atoms_given_template(templates):
    """Select just some atoms in the given template (H2O)."""
    system = templates
    templates = system['template']
    atoms = system['templateatom']
    atoms.current_template = 3  # Acetic acid
    names = []
    x = []
    for row in atoms.atoms('atno', '=', 8, template=2):
        names.append(row['name'])
        x.append(row['x'])
    assert names == ['O']
    assert x == [0.0]


def test_atom_ids(templates):
    """Test getting the atom ids for the current template."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    atoms = system['templateatom']
    ids = atoms.atom_ids()
    assert ids == [*range(4, 12)]


def test_coordinate_system(templates):
    """Test getting and settign the coordinate system."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    atoms = system['templateatom']
    assert atoms.coordinate_system == 'Cartesian'

    with pytest.raises(RuntimeError) as e:
        atoms.coordinate_system = 'fractional'
    assert str(e.value) == 'Templates can only use Cartesian coordinates.'


def test_append_just_coordinates(templates):
    """Test appending just the coordinates."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    atoms = system['templateatom']
    atoms.append(x=[1.0, 2.0, 3.0], y=9.0)
    assert atoms.n_atoms == 11


def test_template_n_bonds(templates):
    """Get the number of bonds in acetic acid template."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    bonds = system['templatebond']
    assert bonds.n_bonds() == 7


def test_template_n_bonds_given_template(templates):
    """Get the number of bonds in water template."""
    system = templates
    templates = system['template']
    templates.set_current_template('acetic acid', type_='molecule')
    bonds = system['templatebond']
    assert bonds.n_bonds(template=2) == 2


def test_append_bonds_error(templates):
    """Append a bond with an invalid atom id type"""
    answer = "'i=a' and 'j=2', the atom indices, must be integers"
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    with pytest.raises(TypeError) as e:
        bonds.append(i=['a'], j=[2])
    assert str(e.value) == answer


def test_append_bonds_scalar(templates):
    """Get the number of bonds in water template."""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    bonds.append(i=1, j=2)
    assert bonds.n_bonds() == 3


def test_append_bonds_invalid_atom_error(templates):
    """Append a bond with an invalid atom id"""
    answer = "Atom i (5) is not in the template."
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    with pytest.raises(ValueError) as e:
        bonds.append(i=5, j=[2])
    assert str(e.value) == answer


def test_append_bonds_invalid_atom_j_error(templates):
    """Append a bond with an invalid atom id"""
    answer = "Atom j (5) is not in the template."
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    with pytest.raises(ValueError) as e:
        bonds.append(i=1, j=5)
    assert str(e.value) == answer


def test_append_bonds_missing_j_error(templates):
    """Append a bond missing the second atom"""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    with pytest.raises(KeyError) as e:
        bonds.append(i=1)
    assert str(e.value) == "'The atoms i & j are required!'"


def test_append_bonds_scalar_i(templates):
    """Append a bond with one i atom and several j's"""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    bonds.append(i=[1], j=[2, 3])
    assert bonds.n_bonds() == 4


def test_append_bonds_scalar_j(templates):
    """Append a bond with one j atom and several i's"""
    system = templates
    templates = system['template']
    templates.set_current_template('H2O', type_='molecule')
    bonds = system['templatebond']
    bonds.append(i=[2, 3], j=[1])
    assert bonds.n_bonds() == 4

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the system class."""

import pytest  # noqa: F401

from molsystem import System

# yapf: disable
#            C       H        H        H        C        =O      O        H
x      = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
y      = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
z      = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
atno   = [6, 1, 1, 1, 6, 8, 8, 1]  # noqa: E221

#         C  H  H  H  C =O  O  H
i_atom = [0, 0, 0, 0, 4, 4, 6]
j_atom = [1, 2, 3, 4, 5, 6, 7]
order  = [1, 1, 1, 1, 1, 2, 1]  # noqa: E221
# yapf: enable


@pytest.fixture()
def bonds():
    """An empty system object
    """
    system = System()
    return system


@pytest.fixture()
def AceticAcid():
    """An system object for an acetic acid molecule
    """
    system = System()
    system['atoms'].append(x=x, y=y, z=z, atno=atno)
    system['bonds'].append(i=i_atom, j=j_atom, order=order)
    return system


def test_construction():
    """Simplest test that we can make a System object"""
    system = System()
    assert str(type(system)) == "<class 'molsystem.system.System'>"


def test_version_empty():
    """Simplest test that we can make a System object"""
    system = System()
    assert system.version == 0


# def test_version_changed():
#     """Simplest test that we can make a System object"""
#     system = molsystem.System()
#     with system as sys:
#         sys.periodicity = 3
#     assert system.version == 1
#
#
# def test_version_unchanged():
#     """Simplest test that we can make a System object"""
#     system = molsystem.System()
#     with system as sys:  # noqa: F841
#         pass
#     assert system.version == 0


def test_atoms_element(AceticAcid):
    assert AceticAcid['atoms'].n_atoms == 8


def test_atoms_attribute(AceticAcid):
    assert AceticAcid.atoms.n_atoms == 8


def test_atoms_with(AceticAcid):
    with AceticAcid.atoms as tmp:
        tmp['x'][0] = 1.0
    assert AceticAcid.atoms.version == 1


def test_system_with(AceticAcid):
    with AceticAcid as tmp:
        tmp.atoms['x'][0] = 1.0
    assert AceticAcid.atoms.version == 1
    assert AceticAcid.version == 1
    assert AceticAcid.bonds.version == 0
    assert AceticAcid.cell.version == 0


def test_system_with_bonds(AceticAcid):
    with AceticAcid as tmp:
        tmp.bonds['order'][0] = 3
    assert AceticAcid.version == 1
    assert AceticAcid.atoms.version == 0
    assert AceticAcid.bonds.version == 1
    assert AceticAcid.cell.version == 0


def test_system_copy_constructor(AceticAcid):
    system = System(AceticAcid)
    assert system.n_atoms == 8
    assert system.n_bonds == 7
    assert system == AceticAcid


def test_str(AceticAcid):
    result = """\
{'atoms':           x       y       z  atno
uid                              
0    1.0797  0.0181 -0.0184   6.0
1    0.5782  3.1376  0.2813   1.0
2    0.7209 -0.6736 -0.7859   1.0
3    0.7052 -0.3143  0.9529   1.0
4    0.5713  1.3899 -0.3161   6.0
5   -0.1323  1.7142 -1.2568   8.0
6    0.9757  2.2970  0.5919   8.0
7    2.1724  0.0161 -0.0306   1.0,
 'bonds':        i    j  order
uid                 
0    0.0  1.0    1.0
1    0.0  2.0    1.0
2    0.0  3.0    1.0
3    0.0  4.0    1.0
4    4.0  5.0    1.0
5    4.0  6.0    2.0
6    6.0  7.0    1.0,
 'cell': {'periodicity': 0, 'cell': [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]}}"""  # noqa: W291, E501

    assert str(AceticAcid) == result


def test_repr(AceticAcid):
    result = """\
{'atoms':           x       y       z  atno
uid                              
0    1.0797  0.0181 -0.0184   6.0
1    0.5782  3.1376  0.2813   1.0
2    0.7209 -0.6736 -0.7859   1.0
3    0.7052 -0.3143  0.9529   1.0
4    0.5713  1.3899 -0.3161   6.0
5   -0.1323  1.7142 -1.2568   8.0
6    0.9757  2.2970  0.5919   8.0
7    2.1724  0.0161 -0.0306   1.0, 'bonds':        i    j  order
uid                 
0    0.0  1.0    1.0
1    0.0  2.0    1.0
2    0.0  3.0    1.0
3    0.0  4.0    1.0
4    4.0  5.0    1.0
5    4.0  6.0    2.0
6    6.0  7.0    1.0, 'cell': {'periodicity': 0, 'cell': [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]}}"""  # noqa: W291, E501

    assert repr(AceticAcid) == result


def test_contains(AceticAcid):
    assert 'atoms' in AceticAcid
    assert 'junk' not in AceticAcid


def test_set_get_del(AceticAcid):
    with AceticAcid as tmp:
        tmp['extra'] = 'extra stuff'
    assert AceticAcid.version == 1
    assert AceticAcid['extra'] == 'extra stuff'
    assert 'extra' in AceticAcid
    with AceticAcid as tmp:
        AceticAcid['extra'] = 'more extra stuff'
    assert AceticAcid.version == 2
    with AceticAcid as tmp:
        del tmp['extra']
    assert AceticAcid.version == 3
    assert 'extra' not in AceticAcid


def test_periodicity(AceticAcid):
    with AceticAcid as tmp:
        tmp.periodicity = 3
    assert AceticAcid.version == 1
    assert AceticAcid.periodicity == 3


def test_no_change(AceticAcid):
    with AceticAcid as tmp:
        tmp.periodicity = 3
        tmp.periodicity = 0


def test_coordinate_type(AceticAcid):
    with AceticAcid as tmp:
        tmp.periodicity = 3
        tmp.coordinate_type = 'frac'
    assert AceticAcid.version == 1
    assert AceticAcid.coordinate_type == 'fractional'

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fixtures for testing the 'molsystem' package."""
import math

import pytest

from molsystem.systems import Systems


def pytest_addoption(parser):
    parser.addoption(
        "--run-timing",
        action="store_true",
        default=False,
        help="run timing tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "timing: mark test as timing to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-timing"):
        # --run-timing given in cli: do not skip timing tests
        return
    skip_timing = pytest.mark.skip(reason="need --run-timing option to run")
    for item in items:
        if "timing" in item.keywords:
            item.add_marker(skip_timing)


def mk_table(system, name='table1'):
    """Create a table with some attributes."""
    table = system.create_table(name)
    table.add_attribute('atno', coltype='int', default=-1)
    for column in ('x', 'y', 'z'):
        table.add_attribute(column, coltype='float')
    return table


@pytest.fixture()
def atoms(system):
    """An empty atoms table."""
    systems = Systems()
    system = systems.create_system('seamm', temporary=True)

    yield system.atoms

    try:
        del systems['seamm']
    except:  # noqa: E722
        print('Caught error deleting the database')
        pass


@pytest.fixture()
def system():
    systems = Systems()
    system = systems.create_system('seamm', temporary=True)

    yield system

    try:
        del systems['seamm']
    except:  # noqa: E722
        print('Caught error deleting the database')
        pass


@pytest.fixture()
def two_systems():
    systems = Systems()

    system1 = systems.create_system('seamm1', temporary=True)
    system2 = systems.create_system('seamm2', temporary=True)

    yield system1, system2

    del systems['seamm1']
    del systems['seamm2']


@pytest.fixture()
def simple_table(system):
    return mk_table(system)


@pytest.fixture()
def two_tables(system):
    return mk_table(system, 'table1'), mk_table(system, 'table2')


@pytest.fixture()
def system_with_two_tables(system):
    mk_table(system, 'table1')
    mk_table(system, 'table2')
    return system


@pytest.fixture()
def two_tables_in_two_systems(two_systems):
    system, system2 = two_systems
    return mk_table(system, 'table1'), mk_table(system2, 'table2')


@pytest.fixture()
def bonds(system):
    """An empty system object
    """
    return system


@pytest.fixture()
def templates(system):
    """A system with a template for water."""
    templates = system['template']
    atoms = system['templateatom']
    bonds = system['templatebond']

    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    name = ['O', 'H1', 'H2']
    i_atom = [0, 0]
    j_atom = [1, 2]

    tid = templates.append(name='H2O', type='molecule')[0]
    templates.current_template = tid

    ids = atoms.append(atno=atno, x=X, y=Y, z=Z, name=name)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    bonds.append(i=j, j=i)  # flipped on purpose so code orders.

    # Acetic acid

    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    X = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    Y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    Z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    symbol = ['C', 'H', 'H', 'H', 'C', 'O', 'O', 'H']
    name = ['C1', 'H1', 'H2', 'H3', 'C', 'O', 'OH', 'H']

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    tid = templates.append(name='acetic acid', type='molecule')[0]
    templates.current_template = tid

    ids = atoms.append(symbol=symbol, x=X, y=Y, z=Z, name=name)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    bonds.append(i=i, j=j, bondorder=order)

    return system


@pytest.fixture()
def AceticAcid(system):
    """An system object for an acetic acid molecule
    """
    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    x = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    atno = [6, 1, 1, 1, 6, 8, 8, 1]  # noqa: E221

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    system.name = 'acetic acid'
    ids = system['atoms'].append(x=x, y=y, z=z, atno=atno)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    system['bonds'].append(i=i, j=j, bondorder=order)

    return system


@pytest.fixture
def disordered(AceticAcid):
    """Two acetic acid molecules with different atom order"""
    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    x = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    atno = [6, 1, 1, 1, 6, 8, 8, 1]  # noqa: E221

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    system = AceticAcid
    ids = system['atoms'].append(
        x=[*reversed(x)],
        y=[*reversed(y)],
        z=[*reversed(z)],
        atno=[*reversed(atno)]
    )

    i = [ids[7 - x] for x in i_atom]
    j = [ids[7 - x] for x in j_atom]

    system['bonds'].append(i=i, j=j, bondorder=order)

    return system


@pytest.fixture()
def vanadium(system):
    """BCC vanadium crystal, without symmetry."""
    system.name = 'BCC Vanadium'
    system.periodicity = 3
    system.coordinate_system = 'fractional'
    system.cell.set_cell(3.03, 3.03, 3.03, 90, 90, 90)
    system.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol='V')
    return system


@pytest.fixture()
def copper(system):
    """FCC copper crystal, without symmetry."""
    x = [0.0, 0.5, 0.5, 0.0]
    y = [0.0, 0.5, 0.0, 0.5]
    z = [0.0, 0.0, 0.5, 0.5]
    system.name = 'FCC Copper'
    system.periodicity = 3
    system.coordinate_system = 'fractional'
    system.cell.set_cell(3.61491, 3.61491, 3.61491, 90, 90, 90)
    system.atoms.append(x=x, y=y, z=z, symbol=['Cu'])
    return system


@pytest.fixture()
def CH3COOH_3H2O(AceticAcid):
    """System with acetic acid and 3 water molecules"""
    system = AceticAcid

    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    i_atom = [0, 0]
    j_atom = [1, 2]

    system = AceticAcid
    system.name = 'acetic acid with 3 waters'

    for no in range(1, 4):
        ids = system['atoms'].append(x=X, y=no * 5.0, z=Z, atno=atno)

        i = [ids[x] for x in i_atom]
        j = [ids[x] for x in j_atom]

        system['bonds'].append(i=i, j=j)

    return system

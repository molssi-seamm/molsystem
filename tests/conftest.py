#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fixtures for testing the 'molsystem' package."""
import os
import shutil
import tempfile

import pytest
from molsystem.system import System


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
def system(tmp_path):
    filepath = tmp_path / 'seamm.db'
    system = System(filename=filepath)
    return system


@pytest.fixture()
def two_systems(tmp_path):
    filepath = tmp_path / 'seamm1.db'
    system1 = System(filename=filepath)

    filepath = tmp_path / 'seamm2.db'
    system2 = System(filename=filepath)

    return system1, system2


@pytest.fixture()
def system2():
    newpath = tempfile.mkdtemp()
    filepath = os.path.join(newpath, 'seamm2.db')

    system = System(filename=filepath)

    yield system

    shutil.rmtree(newpath)


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
    order =  [1, 1, 1, 1, 1, 2, 1]  # noqa: E222
    # yapf: enable

    system['atoms'].append(x=x, y=y, z=z, atno=atno)
    system['bonds'].append(i=i_atom, j=j_atom, order=order)

    return system


@pytest.fixture()
def atoms(system):
    """An empty atoms table."""
    return system['atoms']

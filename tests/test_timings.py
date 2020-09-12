#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest  # noqa: F401
import os
import os.path
import shutil
import tempfile
import time

from molsystem import Systems  # noqa: F401
"""Tests for the System classes."""

natoms = 1000000


@pytest.fixture(scope="module")
def msystem():
    systems = Systems()

    newpath = tempfile.mkdtemp()
    filepath = os.path.join(newpath, 'seamm.db')

    system = systems.create_system('seamm', filename=filepath)

    yield system

    shutil.rmtree(newpath)


@pytest.fixture(scope="module")
def matoms(msystem):
    """An empty atoms table."""
    return msystem['atoms']


@pytest.mark.timing
def test_creation(matoms):
    """Create a system with a million atoms"""
    print(f'\nCreating a system with {natoms} atoms')
    t0 = time.perf_counter()
    rng = numpy.random.default_rng()
    atno = rng.integers(1, 101, size=natoms)
    x = rng.uniform(low=0, high=100, size=natoms)
    y = rng.uniform(low=0, high=100, size=natoms)
    z = rng.uniform(low=0, high=100, size=natoms)
    t1 = time.perf_counter()
    print(f'      setup took {t1-t0:.3} s')

    t0 = time.perf_counter()
    matoms.append(atno=atno.tolist(), x=x.tolist(), y=y.tolist(), z=z.tolist())
    t1 = time.perf_counter()

    print(f'  appending took {t1-t0:.3} s')
    filename = matoms.system.filename
    size = os.path.getsize(filename) / 1024 / 1024
    print(f'       database is {size:.3} MB')

    assert matoms.n_atoms() == natoms


@pytest.mark.timing
def test_access(matoms):
    """Check that we can access the created system"""
    assert matoms.n_atoms() == natoms


@pytest.mark.timing
def test_loop_over_atoms(matoms):
    """Loop over all the atoms (rows)."""
    n = 0
    for atom in matoms.atoms():
        n += 1
    assert n == matoms.n_atoms()


@pytest.mark.timing
def test_selection(matoms):
    """Test selecting a subset of atoms."""
    n = 0
    for atom in matoms.atoms('atno', '==', 6):
        n += 1
    print(f'\nFound {n} carbon atoms in {natoms} atoms')
    print(f'There should be about {int(natoms / 100)}')
    assert atom['atno'] == 6

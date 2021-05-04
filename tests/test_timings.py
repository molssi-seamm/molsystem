#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest  # noqa: F401
import os
import os.path
import time

from molsystem import SystemDB  # noqa: F401

"""Performance tests for molsystem."""

natoms = 1000000


@pytest.fixture(scope="module")
def db_file(tmp_path_factory):
    """Create and temporary directory, and clean it up at the end."""
    dirpath = tmp_path_factory.mktemp("data")
    path = dirpath / "seamm.db"
    db = SystemDB(filename=str(path))
    db.close()

    print(f"temp file: {path}")

    yield str(path)

    path.unlink()
    dirpath.rmdir()


@pytest.fixture()
def disk_db(db_file):
    """Create and return a SystemDB object with the database file."""
    db = SystemDB(filename=db_file)

    yield db

    db.close()


@pytest.fixture(scope="module")
def mem_db():
    """Create and return a SystemDB using an in-memory DB"""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")

    yield db

    db.close()


@pytest.mark.timing
def test_creation(disk_db):
    """Create a system with a million atoms"""
    print(f"\nCreating a system with {natoms} atoms")
    t0 = time.perf_counter()
    system = disk_db.create_system("default")
    configuration = system.create_configuration("default")
    atoms = configuration.atoms
    rng = numpy.random.default_rng()
    atno = rng.integers(1, 101, size=natoms)
    x = rng.uniform(low=0, high=100, size=natoms)
    y = rng.uniform(low=0, high=100, size=natoms)
    z = rng.uniform(low=0, high=100, size=natoms)
    t1 = time.perf_counter()
    print(f"      setup took {t1-t0:.3} s")

    t0 = time.perf_counter()
    atoms.append(atno=atno.tolist(), x=x.tolist(), y=y.tolist(), z=z.tolist())
    t1 = time.perf_counter()

    print(f"  appending took {t1-t0:.3} s")
    filename = disk_db.filename
    size = os.path.getsize(filename) / 1024 / 1024
    print(f"       database is {size:.3} MB")

    assert atoms.n_atoms == natoms


@pytest.mark.timing
def test_access(disk_db):
    """Check that we can access the created system"""
    t0 = time.perf_counter()
    atoms = disk_db.system.configuration.atoms
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    assert atoms.n_atoms == natoms


@pytest.mark.timing
def test_loop_over_atoms(disk_db):
    """Loop over all the atoms (rows)."""
    t0 = time.perf_counter()
    atoms = disk_db.system.configuration.atoms
    n = 0
    for atom in atoms.atoms():
        n += 1
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    assert n == atoms.n_atoms


@pytest.mark.timing
def test_selection(disk_db):
    """Test selecting a subset of atoms."""
    t0 = time.perf_counter()
    atoms = disk_db.system.configuration.atoms
    n = 0
    for atom in atoms.atoms("atno", "==", 6):
        n += 1
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    print(f"\nFound {n} carbon atoms in {natoms} atoms")
    print(f"There should be about {int(natoms / 100)}")
    assert atom["atno"] == 6


@pytest.mark.timing
def test_selection_count(disk_db):
    """Test counting a subset of atoms."""
    t0 = time.perf_counter()
    atoms = disk_db.system.configuration.atoms
    n = atoms.get_n_atoms("atno", "==", 6)
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    print(f"\nFound {n} carbon atoms in {natoms} atoms")
    print(f"There should be about {int(natoms / 100)}")


@pytest.mark.timing
def test_creation_mem(mem_db):
    """Create a system with a million atoms"""
    print(f"\nCreating a system with {natoms} atoms")
    t0 = time.perf_counter()
    system = mem_db.create_system("default")
    configuration = system.create_configuration("default")
    atoms = configuration.atoms
    rng = numpy.random.default_rng()
    atno = rng.integers(1, 101, size=natoms)
    x = rng.uniform(low=0, high=100, size=natoms)
    y = rng.uniform(low=0, high=100, size=natoms)
    z = rng.uniform(low=0, high=100, size=natoms)
    t1 = time.perf_counter()
    print(f"      setup took {t1-t0:.3} s")

    t0 = time.perf_counter()
    atoms.append(atno=atno.tolist(), x=x.tolist(), y=y.tolist(), z=z.tolist())
    t1 = time.perf_counter()

    print(f"  appending took {t1-t0:.3} s")

    assert atoms.n_atoms == natoms


@pytest.mark.timing
def test_access_mem(mem_db):
    """Check that we can access the created system"""
    t0 = time.perf_counter()
    atoms = mem_db.system.configuration.atoms
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    assert atoms.n_atoms == natoms


@pytest.mark.timing
def test_loop_over_atoms_mem(mem_db):
    """Loop over all the atoms (rows)."""
    t0 = time.perf_counter()
    atoms = mem_db.system.configuration.atoms
    n = 0
    for atom in atoms.atoms():
        n += 1
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    assert n == atoms.n_atoms


@pytest.mark.timing
def test_selection_mem(mem_db):
    """Test selecting a subset of atoms."""
    t0 = time.perf_counter()
    atoms = mem_db.system.configuration.atoms
    n = 0
    for atom in atoms.atoms("atno", "==", 6):
        n += 1
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    print(f"\nFound {n} carbon atoms in {natoms} atoms")
    print(f"There should be about {int(natoms / 100)}")
    assert atom["atno"] == 6


@pytest.mark.timing
def test_selection_count_mem(mem_db):
    """Test counting a subset of atoms."""
    t0 = time.perf_counter()
    atoms = mem_db.system.configuration.atoms
    n = atoms.get_n_atoms("atno", "==", 6)
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.3} s")
    print(f"\nFound {n} carbon atoms in {natoms} atoms")
    print(f"There should be about {int(natoms / 100)}")

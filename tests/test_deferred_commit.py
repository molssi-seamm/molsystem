#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Deferred commits: a flowchart step as one transaction (JobConnection)."""

import sqlite3
import subprocess
import sys
import textwrap

import pytest

from molsystem import JobConnection, SystemDB


def count(path, sql):
    """Count with a second, read-only connection: what another process sees."""
    other = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return other.execute(sql).fetchone()[0]
    finally:
        other.close()


@pytest.fixture()
def deferred(tmp_path):
    """A deferring job database on disk with one empty configuration."""
    path = tmp_path / "seamm.db"
    db = SystemDB(filename=f"file:{path}", deferred_commit=True)
    system = db.create_system(name="default")
    system.create_configuration(name="default")
    db.commit_transaction()
    return db, path


def test_connection_class(tmp_path):
    """The job database always uses JobConnection, not deferring by default."""
    db = SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    assert isinstance(db.db, JobConnection)
    assert not db.db.deferring
    assert not db.deferred_commit


def test_not_deferring_commits(tmp_path):
    """Without deferral molsystem's own commits are visible at once."""
    path = tmp_path / "seamm.db"
    db = SystemDB(filename=f"file:{path}")
    db.create_system(name="one")
    assert count(path, "SELECT COUNT(*) FROM system") == 1


def test_deferred_until_commit_transaction(deferred):
    """Writes are invisible to other connections until commit_transaction."""
    db, path = deferred
    configuration = db.system.configuration
    configuration.atoms.append(x=[0.0, 1.0], y=[0.0, 0.0], z=[0.0, 0.0], symbol="H")
    assert configuration.n_atoms == 2
    assert count(path, "SELECT COUNT(*) FROM atom") == 0
    db.commit_transaction()
    assert count(path, "SELECT COUNT(*) FROM atom") == 2


def test_rollback_transaction(deferred):
    """rollback_transaction abandons the step's writes."""
    db, path = deferred
    db.create_system(name="second")
    db.rollback_transaction()
    assert count(path, "SELECT COUNT(*) FROM system") == 1
    assert db.n_systems == 1


def test_killed_process_loses_only_the_open_step(tmp_path):
    """A process killed mid-step keeps earlier steps, loses the running one."""
    path = tmp_path / "seamm.db"
    script = textwrap.dedent(f"""
        import os
        from molsystem import SystemDB

        db = SystemDB(filename="file:{path}", deferred_commit=True)
        system = db.create_system(name="first step")
        system.create_configuration(name="c1")
        db.commit_transaction()

        system = db.create_system(name="second step")
        configuration = system.create_configuration(name="c2")
        configuration.atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol="H")
        os._exit(1)
        """)
    subprocess.run([sys.executable, "-c", script], check=False)
    assert count(path, "SELECT COUNT(*) FROM system") == 1
    db = SystemDB(filename=f"file:{path}")
    assert [s.name for s in db.systems] == ["first step"]


def test_with_configuration_error_restores(deferred):
    """A failing 'with configuration' block restores it without cascades."""
    db, path = deferred
    configuration = db.system.configuration
    configuration.atoms.append(
        x=[0.0, 0.0, 1.1], y=[0.0, 1.1, 0.0], z=[0.0, 0.0, 0.0], symbol=["O", "H", "H"]
    )
    db.commit_transaction()
    n_atoms = configuration.n_atoms
    n_coordinates = count(path, "SELECT COUNT(*) FROM coordinates")

    with pytest.raises(RuntimeError):
        with configuration as tmp:
            tmp.atoms.append(x=[9.0], y=[9.0], z=[9.0], symbol=["He"])
            tmp.periodicity = 3
            raise RuntimeError("fail inside the block")

    # (The configuration object caches its periodicity, with or without
    # deferral, so check the database.)
    assert configuration.n_atoms == n_atoms
    db.commit_transaction()
    assert count(path, "SELECT periodicity FROM configuration") == 0
    assert count(path, "SELECT COUNT(*) FROM atom") == n_atoms
    assert count(path, "SELECT COUNT(*) FROM coordinates") == n_coordinates


def test_with_configuration_success_keeps(deferred):
    """A successful 'with configuration' block keeps its changes, uncommitted."""
    db, path = deferred
    configuration = db.system.configuration
    with configuration as tmp:
        tmp.periodicity = 3
        tmp.coordinate_system = "fractional"
        tmp.cell.parameters = (3.03, 3.03, 3.03, 90, 90, 90)
        tmp.atoms.append(x=[0.0, 0.5], y=[0.0, 0.5], z=[0.0, 0.5], symbol="V")
    assert configuration.atoms.n_atoms == 2
    assert configuration.version == 1
    assert count(path, "SELECT COUNT(*) FROM atom") == 0
    db.commit_transaction()
    assert count(path, "SELECT COUNT(*) FROM atom") == 2


def test_nested_with_blocks(deferred):
    """An inner failing block rolls back only itself."""
    db, path = deferred
    atoms = db.system.configuration.atoms
    with atoms:
        atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol="H")
        with pytest.raises(ValueError):
            with atoms:
                atoms.append(x=[1.0], y=[0.0], z=[0.0], symbol="H")
                raise ValueError("inner")
    assert atoms.n_atoms == 1


def test_delete_column_deferred(deferred):
    """Dropping a column does not commit the deferred transaction."""
    db, path = deferred
    atoms = db.system.configuration.atoms
    atoms.add_attribute("charge", coltype="float", default=0.0)
    atoms.append(x=[0.0], y=[0.0], z=[0.0], symbol="H")
    assert "charge" in db["atom"].attributes
    del db["atom"]["charge"]
    assert "charge" not in db["atom"].attributes
    assert count(path, "SELECT COUNT(*) FROM atom") == 0
    db.commit_transaction()
    assert count(path, "SELECT COUNT(*) FROM atom") == 1


def test_first_write_of_a_step_is_create_table(tmp_path):
    """DDL right after a commit is still part of the step's transaction."""
    path = tmp_path / "seamm.db"
    script = textwrap.dedent(f"""
        import os
        from molsystem import SystemDB

        db = SystemDB(filename="file:{path}", deferred_commit=True)
        db.create_system(name="first step")
        db.commit_transaction()

        db.db.execute("CREATE TABLE made_by_the_second_step (x)")
        db.user_tables.create("table of the second step", columns=[("x", "float")])
        os._exit(1)
        """)
    subprocess.run([sys.executable, "-c", script], check=False)
    other = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    names = [r[0] for r in other.execute("SELECT name FROM sqlite_master")]
    other.close()
    assert "made_by_the_second_step" not in names
    assert not any(name.startswith("table_") for name in names)
    assert count(path, "SELECT COUNT(*) FROM system") == 1


def test_always_in_a_transaction_while_deferring(deferred):
    db, path = deferred
    assert db.db.in_transaction
    db.commit_transaction()
    assert db.db.in_transaction
    db.rollback_transaction()
    assert db.db.in_transaction
    db.deferred_commit = False
    db.commit_transaction()
    assert not db.db.in_transaction

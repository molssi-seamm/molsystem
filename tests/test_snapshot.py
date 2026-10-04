#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Snapshots for parallel loop iterations, and their merge (phase 6)."""

import sqlite3

import pytest

from molsystem import SystemDB
from molsystem.snapshot import MergeConflict, baseline, merge, snapshot

WATER = dict(
    x=[0.0, 0.0, 0.76], y=[0.0, 0.76, 0.0], z=[0.0, 0.0, 0.0], symbol=["O", "H", "H"]
)


def make_job(path):
    """A job database: three systems of one configuration, a property, a table."""
    db = SystemDB(filename=f"file:{path}", deferred_commit=True)
    db.properties.add("energy", units="kJ/mol")
    for name in ("first", "second", "third"):
        system = db.create_system(name=name)
        configuration = system.create_configuration(name=name)
        configuration.atoms.append(**WATER)
        configuration.bonds.append(i=[1, 1], j=[2, 3])
        db.properties.put(configuration.id, "energy", -1.0)
    table = db.user_tables.create(
        "results",
        columns=[("name", "string", ""), ("E", "float", None)],
        index_column="name",
    )
    table.append_rows([{"name": "first"}, {"name": "second"}, {"name": "third"}])
    db.commit_transaction()
    return db


def configuration_named(db, name):
    return db.get_system(name).configuration


def iteration(tmp_path, job, name, k):
    """Snapshot for the configuration of system ``name``; return the child db."""
    cid = configuration_named(job, name).id
    child_path = tmp_path / f"iter_{k}" / "seamm.db"
    snapshot(tmp_path / "seamm.db", child_path, configurations=[cid])
    baseline(child_path, tmp_path / f"iter_{k}" / "baseline.db")
    return SystemDB(filename=f"file:{child_path}"), child_path


def committed(path, sql):
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return db.execute(sql).fetchall()
    finally:
        db.close()


def test_snapshot_holds_only_the_selected_structures(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "second", 1)
    assert [s.name for s in child.systems] == ["second"]
    configuration = child.system.configuration
    assert configuration.id == configuration_named(job, "second").id  # same ids
    assert configuration.n_atoms == 3
    assert configuration.bonds.n_bonds == 2
    assert child.properties.get(configuration.id, "energy")["energy"]["value"] == -1.0
    rows = child.user_tables["results"].row_ids()
    assert len(rows) == 3  # all tables
    assert committed(path, "SELECT COUNT(*) FROM _table_changes") == [(0,)]


def test_snapshot_reads_only_committed_state(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    job.create_system(name="uncommitted").create_configuration(name="x")
    cid = configuration_named(job, "first").id
    snapshot(tmp_path / "seamm.db", tmp_path / "c.db", configurations=[cid], whole=True)
    names = [r[0] for r in committed(tmp_path / "c.db", "SELECT name FROM system")]
    assert "uncommitted" not in names
    assert len(names) == 3


def test_merge_brings_back_structures_properties_and_rows(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "second", 1)
    configuration = child.system.configuration
    cid = configuration.id
    # What a body might do
    configuration.atoms.set_coordinates([[0, 0, 0], [0, 0.8, 0], [0.8, 0, 0]])
    child.properties.add("dipole", units="debye")
    child.properties.put(cid, "dipole", 1.85)
    child.properties.put(cid, "energy", -2.0)
    new = child.create_system(name="made in the loop")
    made = new.create_configuration(name="made")
    made.atoms.append(x=[0.0, 0.74], y=[0.0, 0.0], z=[0.0, 0.0], symbol=["H", "H"])
    child.properties.put(made.id, "energy", -0.5)
    table = child.user_tables["results"]
    table.set_cell(table.rowid_at(1), "E", -2.0)  # the row of "second"
    table.add_column("dipole", "float", None)
    table.set_cell(table.rowid_at(1), "dipole", 1.85)
    table.append_row(name="made in the loop", E=-0.5)
    child.db.commit()
    child.close()

    result = merge(job, path, tmp_path / "iter_1" / "baseline.db", iteration=1)
    job.commit_transaction()

    second = configuration_named(job, "second")
    assert second.atoms.get_coordinates()[2] == pytest.approx([0.8, 0.0, 0.0])
    assert job.properties.get(second.id, "dipole")["dipole"]["value"] == 1.85
    assert job.properties.get(second.id, "energy")["energy"]["value"] == -2.0
    names = [s.name for s in job.systems]
    assert names == ["first", "second", "third", "made in the loop"]
    made_in_job = job.get_system("made in the loop").configuration
    assert made_in_job.n_atoms == 2
    assert made_in_job.atoms.symbols == ["H", "H"]
    assert job.properties.get(made_in_job.id, "energy")["energy"]["value"] == -0.5
    results = job.user_tables["results"]
    assert results.columns == ["name", "E", "dipole"]
    assert results.get_cell(results.rowid_at(1), "E") == -2.0
    assert results.get_cell(results.rowid_at(3), "name") == "made in the loop"
    # The untouched structures are as they were
    first = configuration_named(job, "first")
    assert first.atoms.get_coordinates()[1] == pytest.approx([0.0, 0.76, 0.0])
    assert result["maps"]["system"]


def test_two_iterations_new_structures_and_rows_do_not_collide(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    state = {}
    paths = []
    for k, name in enumerate(("first", "third"), start=1):
        child, path = iteration(tmp_path, job, name, k)
        new = child.create_system(name=f"new {k}")
        c = new.create_configuration(name=f"new {k}")
        c.atoms.append(x=[float(k)], y=[0.0], z=[0.0], symbol=["He"])
        child.user_tables["results"].append_row(name=f"row {k}", E=float(k))
        child.db.commit()
        child.close()
        paths.append(path)
    for k, path in enumerate(paths, start=1):
        merge(job, path, tmp_path / f"iter_{k}" / "baseline.db", state, iteration=k)
    job.commit_transaction()
    assert [s.name for s in job.systems][-2:] == ["new 1", "new 2"]
    assert job.get_system("new 2").configuration.atoms.get_coordinates()[0][0] == 2.0
    results = job.user_tables["results"]
    assert [results.get_cell(r, "name") for r in results.row_ids()][-2:] == [
        "row 1",
        "row 2",
    ]
    assert job.n_systems == 5


def test_same_cell_from_two_iterations(tmp_path):
    for later_wins in (False, True):
        directory = tmp_path / str(later_wins)
        directory.mkdir()
        job = make_job(directory / "seamm.db")
        state = {}
        paths = []
        for k, name in enumerate(("first", "second"), start=1):
            child, path = iteration(directory, job, name, k)
            table = child.user_tables["results"]
            table.set_cell(table.rowid_at(2), "E", float(k))  # both: "third"'s row
            child.db.commit()
            child.close()
            paths.append(path)
        merge(job, paths[0], directory / "iter_1" / "baseline.db", state, 1)
        if later_wins:
            merge(
                job,
                paths[1],
                directory / "iter_2" / "baseline.db",
                state,
                2,
                later_wins=True,
            )
            results = job.user_tables["results"]
            assert results.get_cell(results.rowid_at(2), "E") == 2.0
        else:
            with pytest.raises(MergeConflict, match="Iterations 1 and 2"):
                merge(job, paths[1], directory / "iter_2" / "baseline.db", state, 2)
        job.rollback_transaction()


def test_same_index_appended_by_two_iterations(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    state = {}
    paths = []
    for k, name in enumerate(("first", "second"), start=1):
        child, path = iteration(tmp_path, job, name, k)
        child.user_tables["results"].append_row(name="shared", E=float(k))
        child.db.commit()
        child.close()
        paths.append(path)
    merge(job, paths[0], tmp_path / "iter_1" / "baseline.db", state, 1)
    with pytest.raises(MergeConflict, match="shared"):
        merge(job, paths[1], tmp_path / "iter_2" / "baseline.db", state, 2)


def test_current_row_and_new_table(tmp_path):
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "first", 1)
    table = child.user_tables["results"]
    table.current_row = table.rowid_at(0)  # it was the last row
    child.user_tables.create("made", columns=[("x", "float", 0.0)])
    child.user_tables["made"].append_row(x=1.5)
    child.db.commit()
    child.close()
    result = merge(job, path, tmp_path / "iter_1" / "baseline.db", {}, 1)
    job.commit_transaction()
    assert result["current_rows"]["results"] == job.user_tables["results"].rowid_at(0)
    made = job.user_tables["made"]
    assert [made.get_cell(r, "x") for r in made.row_ids()] == [1.5]


def test_merge_is_one_transaction(tmp_path):
    """A merge rolled back leaves the job's database as it was."""
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "first", 1)
    child.create_system(name="new").create_configuration(name="new")
    child.db.commit()
    child.close()
    merge(job, path, tmp_path / "iter_1" / "baseline.db", {}, 1)
    job.rollback_transaction()
    assert [s.name for s in job.systems] == ["first", "second", "third"]


def test_new_configuration_of_an_existing_system(tmp_path):
    """A conformer shares its system's atoms; the merge must not copy them."""
    job = make_job(tmp_path / "seamm.db")
    n_atoms = committed(tmp_path / "seamm.db", "SELECT COUNT(*) FROM atom")[0][0]
    child, path = iteration(tmp_path, job, "second", 1)
    system = child.system
    conformer = system.copy_configuration(name="optimized", make_current=True)
    conformer.atoms.set_coordinates([[0, 0, 0.1], [0, 0.8, 0], [0.8, 0, 0]])
    child.db.commit()
    child.close()
    merge(job, path, tmp_path / "iter_1" / "baseline.db", {}, 1)
    job.commit_transaction()
    second = job.get_system("second")
    assert [c.name for c in second.configurations] == ["second", "optimized"]
    optimized = second.configuration  # it became the default
    assert optimized.name == "optimized"
    assert optimized.atoms.get_coordinates()[0] == pytest.approx([0, 0, 0.1])
    original = second.get_configuration("second")
    assert original.atoms.get_coordinates()[0] == pytest.approx([0, 0, 0])
    assert optimized.atomset == original.atomset
    assert (
        committed(tmp_path / "seamm.db", "SELECT COUNT(*) FROM atom")[0][0] == n_atoms
    )


def test_attribute_added_by_an_iteration(tmp_path):
    """A step adds a per-atom attribute (MOPAC's charges): a new column."""
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "second", 1)
    configuration = child.system.configuration
    configuration.atoms.add_attribute(
        "charge", coltype="float", configuration_dependent=True
    )
    configuration.atoms["charge"] = [-0.8, 0.4, 0.4]
    child.db.commit()
    child.close()
    merge(job, path, tmp_path / "iter_1" / "baseline.db", {}, 1)
    job.commit_transaction()
    second = configuration_named(job, "second")
    assert second.atoms.get_column_data("charge") == pytest.approx([-0.8, 0.4, 0.4])
    first = configuration_named(job, "first")
    assert first.atoms.get_column_data("charge") == [None, None, None]


def test_row_appended_before_its_column_was_added(tmp_path):
    """A row appended, then a step adds a column and fills it (MOPAC results)."""
    job = make_job(tmp_path / "seamm.db")
    child, path = iteration(tmp_path, job, "first", 1)
    table = child.user_tables["results"]
    row = table.append_row(name="new")
    table.add_column("dipole", "float", None)
    table.set_cell(row, "dipole", 1.5)
    child.db.commit()
    child.close()
    merge(job, path, tmp_path / "iter_1" / "baseline.db", {}, 1)
    job.commit_transaction()
    results = job.user_tables["results"]
    assert results.get_cell(results.rowid_at(3), "dipole") == 1.5

# -*- coding: utf-8 -*-

"""Snapshots of a job database for parallel loop iterations, and their merge.

A parallel Loop runs each iteration in its own evaluator with its own database: a
*snapshot* of the job's database holding the configurations the iteration works on
(with their own ids), all the user tables and the property definitions. When the
iteration has finished, :func:`merge` brings what it did back into the job's
database: new and changed structures (new ones get new ids), property values, and
table rows replayed from the iteration's change journal.

Both work on committed state through their own connections, never through the job
database's (possibly deferring) connection for reading the other file.
"""

import logging
from pathlib import Path
import shutil
import sqlite3

logger = logging.getLogger(__name__)

# Tables a snapshot leaves empty: they belong to the run, not to the data.
RUN_TABLES = ("_table_changes", "_checkpoint", "_checkpoint_variables")

# Structure tables keyed by id, in the order their foreign keys allow.
ID_TABLES = (
    "symmetry",
    "cell",
    "atom",
    "atomset",
    "bond",
    "bondset",
    "system",
    "configuration",
    "template",
    "subset",
)

# Foreign keys of the structure tables: column -> referenced table.
FOREIGN_KEYS = {
    "configuration": {
        "system": "system",
        "symmetry": "symmetry",
        "cell": "cell",
        "atomset": "atomset",
        "bondset": "bondset",
    },
    "system": {"default_configuration": "configuration"},
    "template": {"configuration": "configuration"},
    "bond": {"i": "atom", "j": "atom"},
    "atomset_atom": {"atomset": "atomset", "atom": "atom"},
    "bondset_bond": {"bondset": "bondset", "bond": "bond"},
    "coordinates": {"configuration": "configuration", "atom": "atom"},
    "velocities": {"configuration": "configuration", "atom": "atom"},
    "gradients": {"configuration": "configuration", "atom": "atom"},
    "subset": {"configuration": "configuration", "template": "template"},
    "subset_atom": {"subset": "subset", "atom": "atom", "templateatom": "atom"},
}

# Tables of (configuration, atom) data, and link tables, keyed by these columns.
ROW_TABLES = {
    "atomset_atom": ("atomset", "atom"),
    "bondset_bond": ("bondset", "bond"),
    "coordinates": ("configuration", "atom"),
    "velocities": ("configuration", "atom"),
    "gradients": ("configuration", "atom"),
    "subset_atom": ("subset", "atom"),
}

PROPERTY_DATA = ("float_data", "int_data", "str_data", "json_data")


def _tables(db, schema="main"):
    """The tables of a database: name -> CREATE statement."""
    return {
        name: sql
        for name, sql in db.execute(
            f"SELECT name, sql FROM {schema}.sqlite_master WHERE type = 'table'"
            " AND name NOT LIKE 'sqlite_%'"
        )
    }


def _columns(db, table, schema="main"):
    return [row[1] for row in db.execute(f'PRAGMA {schema}.table_info("{table}")')]


def snapshot(source, target, configurations=(), whole=False):
    """Write a snapshot of a job database for one parallel loop iteration.

    Parameters
    ----------
    source : str or pathlib.Path
        The job database. Only its committed state is read, through a separate
        read-only connection.
    target : str or pathlib.Path
        The new database file (replaced if it exists).
    configurations : [int]
        The ids of the configurations the iteration works on. Their systems,
        atoms, bonds, cells, symmetry, coordinates, velocities, gradients, subsets
        and property values are copied with their own ids; so are the
        configurations of the templates.
    whole : bool
        Copy the whole database instead (the Loop's "whole database" option).

    The change journal and the checkpoint tables are left empty.
    """
    source = Path(source)
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("", "-wal", "-shm"):
        Path(str(target) + suffix).unlink(missing_ok=True)

    if whole:
        src = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
        dst = sqlite3.connect(str(target))
        try:
            src.backup(dst)
            for table in RUN_TABLES:
                if table in _tables(dst):
                    dst.execute(f'DELETE FROM "{table}"')
            dst.commit()
        finally:
            src.close()
            dst.close()
        return

    dst = sqlite3.connect(str(target), uri=True)
    try:
        dst.execute("PRAGMA foreign_keys = OFF")
        dst.execute(f"ATTACH DATABASE 'file:{source}?mode=ro' AS src")
        tables = _tables(dst, "src")
        # The same schema, including columns steps added to atom, coordinates, ...
        for name, sql in tables.items():
            dst.execute(sql)
        for (sql,) in dst.execute(
            "SELECT sql FROM src.sqlite_master WHERE type = 'index' AND sql IS NOT NULL"
        ).fetchall():
            dst.execute(sql)

        # The ids to copy
        dst.execute("CREATE TEMP TABLE sel_conf (id INTEGER PRIMARY KEY)")
        dst.executemany(
            "INSERT OR IGNORE INTO sel_conf VALUES (?)",
            [(int(c),) for c in configurations],
        )
        if "template" in tables:
            dst.execute(
                "INSERT OR IGNORE INTO sel_conf SELECT configuration FROM src.template"
                " WHERE configuration IS NOT NULL"
            )
        selections = {
            "configuration": "SELECT id FROM sel_conf",
            "system": "SELECT system FROM src.configuration WHERE id IN sel_conf",
            "cell": "SELECT cell FROM src.configuration WHERE id IN sel_conf",
            "symmetry": "SELECT symmetry FROM src.configuration WHERE id IN sel_conf",
            "atomset": "SELECT atomset FROM src.configuration WHERE id IN sel_conf",
            "bondset": "SELECT bondset FROM src.configuration WHERE id IN sel_conf",
        }
        for table, select in selections.items():
            if table in tables:
                dst.execute(
                    f'INSERT INTO main."{table}" SELECT * FROM src."{table}"'
                    f" WHERE id IN ({select})"
                )
        if "atomset_atom" in tables:
            dst.execute(
                "INSERT INTO main.atomset_atom SELECT * FROM src.atomset_atom"
                " WHERE atomset IN (SELECT id FROM main.atomset)"
            )
            dst.execute(
                "INSERT INTO main.atom SELECT * FROM src.atom"
                " WHERE id IN (SELECT atom FROM main.atomset_atom)"
            )
        if "bondset_bond" in tables:
            dst.execute(
                "INSERT INTO main.bondset_bond SELECT * FROM src.bondset_bond"
                " WHERE bondset IN (SELECT id FROM main.bondset)"
            )
            dst.execute(
                "INSERT INTO main.bond SELECT * FROM src.bond"
                " WHERE id IN (SELECT bond FROM main.bondset_bond)"
            )
        for table in ("coordinates", "velocities", "gradients", "subset"):
            if table in tables:
                dst.execute(
                    f'INSERT INTO main."{table}" SELECT * FROM src."{table}"'
                    " WHERE configuration IN sel_conf"
                )
        if "subset_atom" in tables:
            dst.execute(
                "INSERT INTO main.subset_atom SELECT * FROM src.subset_atom"
                " WHERE subset IN (SELECT id FROM main.subset)"
            )
        for table in PROPERTY_DATA:
            if table in tables:
                dst.execute(
                    f'INSERT INTO main."{table}" SELECT * FROM src."{table}"'
                    " WHERE configuration IN sel_conf OR (configuration IS NULL AND"
                    " system IN (SELECT id FROM main.system))"
                )

        # Whole tables: shared definitions and the user tables
        done = {
            *selections,
            "atom",
            "atomset_atom",
            "bond",
            "bondset_bond",
            "coordinates",
            "velocities",
            "gradients",
            "subset",
            "subset_atom",
            *PROPERTY_DATA,
            *RUN_TABLES,
        }
        for table in tables:
            if table in done:
                continue
            if table.endswith("_collection_system") or table.endswith(
                "_collection_configuration"
            ):
                continue  # links to structures that may not be here
            dst.execute(f'INSERT INTO main."{table}" SELECT * FROM src."{table}"')
        if dst.execute(
            "SELECT 1 FROM src.sqlite_master WHERE name = 'sqlite_sequence'"
        ).fetchone():
            dst.execute("DELETE FROM main.sqlite_sequence")
            dst.execute(
                "INSERT INTO main.sqlite_sequence SELECT * FROM src.sqlite_sequence"
            )
        dst.commit()
        dst.execute("DETACH DATABASE src")
    finally:
        dst.close()


def baseline(snapshot_path, baseline_path):
    """Keep an untouched copy of a snapshot, read-only, to diff the merge against."""
    baseline_path = Path(baseline_path)
    baseline_path.unlink(missing_ok=True)
    shutil.copyfile(snapshot_path, baseline_path)
    baseline_path.chmod(0o444)


class MergeConflict(RuntimeError):
    """Two iterations of a parallel loop changed the same thing."""


def merge(target, source, baseline_path, state=None, iteration=None, later_wins=False):
    """Bring what a parallel loop iteration did back into the job's database.

    Parameters
    ----------
    target : molsystem.SystemDB
        The job's database. Writes go through its connection, inside its
        transaction; the caller commits (with the checkpoint).
    source : str or pathlib.Path
        The iteration's database, read through its own read-only connection.
    baseline_path : str or pathlib.Path
        The untouched snapshot the iteration started from.
    state : dict
        Carried from one iteration's merge to the next in the same loop: what
        earlier iterations changed, to detect two changing the same thing.
    iteration : any
        The iteration, for messages.
    later_wins : bool
        If two iterations change the same table cell or structure, let the later
        one win (with a warning) instead of raising MergeConflict.

    Returns
    -------
    dict
        ``maps`` (table -> {iteration id: job id}), ``current_rows`` (table name ->
        the job's row the iteration left current, for tables whose current row it
        moved), ``exported`` (tables whose export metadata it changed).
    """
    if state is None:
        state = {}
    touched = state.setdefault("touched", {})
    child = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
    base = sqlite3.connect(f"file:{baseline_path}?mode=ro", uri=True)
    db = target.db
    try:
        maps = {}
        _merge_properties(db, child, base, maps)
        _merge_structures(db, child, base, maps, touched, iteration, later_wins)
        result = _merge_tables(
            target, source, base, maps, touched, iteration, later_wins
        )
    finally:
        child.close()
        base.close()
    result["maps"] = maps
    return result


def _conflict(touched, key, iteration, later_wins, what):
    earlier = touched.get(key)
    if earlier is not None and earlier != iteration:
        message = (
            f"Iterations {earlier} and {iteration} of the parallel loop both changed "
            f"{what}."
        )
        if not later_wins:
            raise MergeConflict(message + " (The Loop can let the later one win.)")
        logger.warning(message + " The later one wins.")
    touched[key] = iteration


def _rows(db, table, columns):
    selected = ", ".join(f'"{c}"' for c in columns)
    return db.execute(f'SELECT {selected} FROM "{table}"').fetchall()


def _merge_properties(db, child, base, maps):
    """Property definitions by name; new ones created in the job's database."""
    if "property" not in _tables(child):
        return
    mapping = {}
    for pid, name, ptype, units, description in child.execute(
        "SELECT id, name, type, units, description FROM property"
    ):
        row = db.execute("SELECT id FROM property WHERE name = ?", (name,)).fetchone()
        if row is None:
            cursor = db.execute(
                "INSERT INTO property (name, type, units, description)"
                " VALUES (?, ?, ?, ?)",
                (name, ptype, units, description),
            )
            mapping[pid] = cursor.lastrowid
        else:
            mapping[pid] = row[0]
    maps["property"] = mapping


def _merge_structures(db, child, base, maps, touched, iteration, later_wins):
    child_tables = _tables(child)
    target_tables = _tables(db)
    deferred = []  # (table, job id, column, iteration id) for forward references

    def mapped(table, column, value, row_table):
        ref = FOREIGN_KEYS.get(row_table, {}).get(column)
        if ref is None or value is None:
            return value
        if ref in maps and value in maps[ref]:
            return maps[ref][value]
        return None  # not mapped yet: filled in later

    # Tables keyed by id
    for table in ID_TABLES:
        if table not in child_tables or table not in target_tables:
            continue
        columns = _columns(child, table)
        new_rows = {row[0]: row for row in _rows(child, table, columns)}
        old_rows = {row[0]: row for row in _rows(base, table, columns)}
        mapping = maps.setdefault(table, {})
        others = columns[1:]
        for _id, row in new_rows.items():
            if _id in old_rows:
                mapping[_id] = _id
        for _id, row in new_rows.items():
            values = [mapped(table, c, v, table) for c, v in zip(others, row[1:])]
            for c, v, value in zip(others, row[1:], values):
                if value is None and v is not None:
                    deferred.append((table, _id, c, v))
            if _id in old_rows:
                if row == old_rows[_id]:
                    continue
                _conflict(
                    touched, (table, _id), iteration, later_wins, f"{table} {_id}"
                )
                assignments = ", ".join(f'"{c}" = ?' for c in others)
                db.execute(
                    f'UPDATE "{table}" SET {assignments} WHERE id = ?', (*values, _id)
                )
            else:
                if len(others) == 0:  # e.g. atomset, bondset: only an id
                    cursor = db.execute(f'INSERT INTO "{table}" DEFAULT VALUES')
                else:
                    names = ", ".join(f'"{c}"' for c in others)
                    places = ", ".join("?" * len(others))
                    cursor = db.execute(
                        f'INSERT INTO "{table}" ({names}) VALUES ({places})', values
                    )
                mapping[_id] = cursor.lastrowid
        # Structures the iteration deleted
        if table in ("system", "configuration"):
            for _id in old_rows:
                if _id not in new_rows:
                    _conflict(
                        touched, (table, _id), iteration, later_wins, f"{table} {_id}"
                    )
                    db.execute(f'DELETE FROM "{table}" WHERE id = ?', (_id,))

    # References to rows inserted after the row that refers to them
    for table, _id, column, value in deferred:
        ref = FOREIGN_KEYS[table][column]
        job_id = maps[table][_id]
        db.execute(
            f'UPDATE "{table}" SET "{column}" = ? WHERE id = ?',
            (maps.get(ref, {}).get(value), job_id),
        )

    # Link and per-atom tables, keyed by columns
    for table, key in ROW_TABLES.items():
        if table not in child_tables or table not in target_tables:
            continue
        columns = _columns(child, table)
        index = [columns.index(k) for k in key]

        def job_row(row):
            return tuple(mapped(table, c, v, table) for c, v in zip(columns, row))

        new_rows = {
            tuple(row[i] for i in index): row for row in _rows(child, table, columns)
        }
        old_rows = {
            tuple(row[i] for i in index): row for row in _rows(base, table, columns)
        }
        where = " AND ".join(f'"{k}" = ?' for k in key)
        for k, row in new_rows.items():
            if old_rows.get(k) == row:
                continue
            values = job_row(row)
            job_key = tuple(values[i] for i in index)
            db.execute(f'DELETE FROM "{table}" WHERE {where}', job_key)
            names = ", ".join(f'"{c}"' for c in columns)
            places = ", ".join("?" * len(columns))
            db.execute(f'INSERT INTO "{table}" ({names}) VALUES ({places})', values)
        for k, row in old_rows.items():
            if k not in new_rows:
                job_key = tuple(mapped(table, c, v, table) for c, v in zip(key, k))
                db.execute(f'DELETE FROM "{table}" WHERE {where}', job_key)

    # Property values, keyed by (configuration, system, property)
    pmap = maps.get("property", {})
    for table in PROPERTY_DATA:
        if table not in child_tables or table not in target_tables:
            continue
        columns = _columns(child, table)  # id, configuration, system, property, value

        def key_of(row):
            r = dict(zip(columns, row))
            return (r["configuration"], r["system"], r["property"])

        new_rows = {key_of(r): r for r in _rows(child, table, columns)}
        old_rows = {key_of(r): r for r in _rows(base, table, columns)}
        for k, row in new_rows.items():
            old = old_rows.get(k)
            if old is not None and old[1:] == row[1:]:
                continue
            configuration, system, prop = k
            job_conf = maps.get("configuration", {}).get(configuration, configuration)
            job_sys = maps.get("system", {}).get(system, system)
            job_prop = pmap.get(prop, prop)
            value = dict(zip(columns, row))["value"]
            if job_conf is None:
                db.execute(
                    f'DELETE FROM "{table}" WHERE configuration IS NULL AND system = ?'
                    " AND property = ?",
                    (job_sys, job_prop),
                )
            else:
                db.execute(
                    f'DELETE FROM "{table}" WHERE configuration = ? AND property = ?',
                    (job_conf, job_prop),
                )
            db.execute(
                f'INSERT INTO "{table}" (configuration, system, property, value)'
                " VALUES (?, ?, ?, ?)",
                (job_conf, job_sys, job_prop, value),
            )


def _merge_tables(target, source, base, maps, touched, iteration, later_wins):
    """Replay the iteration's table journal into the job's tables."""
    from .system_db import SystemDB  # noqa: F811

    result = {"current_rows": {}, "exported": []}
    child = SystemDB(filename=f"file:{source}?mode=ro")
    try:
        if "_tables" not in _tables(child.db):
            return result
        tables = target.user_tables
        ctables = child.user_tables
        journal = ctables.journal()
        base_rows = {}
        base_registry = {}
        if "_tables" in _tables(base):
            for name, sql_name, current, metadata in base.execute(
                "SELECT name, sql_name, current_row, metadata FROM _tables"
            ):
                base_registry[name] = (current, metadata)
                base_rows[name] = {
                    r[0] for r in base.execute(f'SELECT "__rowid__" FROM "{sql_name}"')
                }
        rowmap = {}
        for entry in journal:
            name, op = entry["table"], entry["op"]
            row, column = entry["row"], entry["column"]
            if op == "create":
                if name not in tables:
                    source_table = ctables[name]
                    tables.create(
                        name,
                        columns=[
                            (c["name"], c["type"], c["default"])
                            for c in source_table.column_definitions
                        ],
                        index_column=source_table.index_column,
                        metadata=source_table.metadata,
                    )
            elif op == "drop":
                if name in tables and name not in ctables:
                    tables.delete(name)
            elif op == "add_column":
                definition = ctables[name]._definition(column)
                table = tables[name]
                if column in table.columns:
                    mine = table._definition(column)
                    if (mine["type"], mine["default"]) != (
                        definition["type"],
                        definition["default"],
                    ):
                        raise MergeConflict(
                            f"Iteration {iteration} gave the column '{column}' of the "
                            f"table '{name}' type {definition['type']} and default "
                            f"{definition['default']!r}, but the job's table has "
                            f"{mine['type']} and {mine['default']!r}."
                        )
                else:
                    table.add_column(column, definition["type"], definition["default"])
            elif op == "append":
                if name not in ctables:
                    continue
                values = ctables[name].get_row(row)
                table = tables[name]
                index = table.index_column
                if index is not None and index in values:
                    existing = table.find(index, values[index])
                    if len(existing) > 0:
                        _conflict(
                            touched,
                            ("table", name, "index", values[index]),
                            iteration,
                            later_wins,
                            f"the row '{values[index]}' of the table '{name}'",
                        )
                        for c, v in values.items():
                            table.set_cell(existing[0], c, v)
                        rowmap[(name, row)] = existing[0]
                        continue
                touched[("table", name, "index", values.get(index))] = iteration
                (rowmap[(name, row)],) = table.append_rows([values], move_current=False)
            elif op == "set":
                if name not in ctables or not ctables[name].has_row(row):
                    continue
                if (name, row) in rowmap:
                    job_row = rowmap[(name, row)]
                elif row in base_rows.get(name, ()):
                    job_row = row
                    _conflict(
                        touched,
                        ("cell", name, row, column),
                        iteration,
                        later_wins,
                        f"the cell '{column}' of row {row} of the table '{name}'",
                    )
                else:
                    continue
                tables[name].set_cell(
                    job_row, column, ctables[name].get_cell(row, column)
                )
        # Current rows the iteration moved, and exports it made
        for name in ctables:
            ctable = ctables[name]
            current = ctable.current_row
            before = base_registry.get(name, (None, None))
            if current is not None and current != before[0]:
                result["current_rows"][name] = rowmap.get((name, current), current)
            metadata = ctable.metadata
            if metadata.get("filename") and (
                before[1] is None or metadata != _loads(before[1])
            ):
                result["exported"].append(name)
    finally:
        child.close()
    return result


def _loads(text):
    import json

    try:
        return json.loads(text)
    except (TypeError, ValueError):
        return {}

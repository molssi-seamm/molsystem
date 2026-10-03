#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the user tables (flowchart results tables) in a SystemDB."""

import math

import numpy as np
import pandas
import pytest

from molsystem import SystemDB


@pytest.fixture()
def tables(empty_db):
    return empty_db.user_tables


def test_empty(tables):
    """A new database has no user tables and no registry."""
    assert tables.names == []
    assert "x" not in tables
    assert "_tables" not in tables.system_db


def test_create(tables):
    table = tables.create(
        "Energies", columns=[("Name", "string", None), ("E (kJ/mol)", "float", None)]
    )
    assert tables.names == ["Energies"]
    assert table.columns == ["Name", "E (kJ/mol)"]
    assert table.sql_name == "table_1"
    assert table.n_rows == 0
    assert table.current_row is None
    assert table.column_type("E (kJ/mol)") == "float"
    assert math.isnan(table.default("E (kJ/mol)"))
    assert table.default("Name") == ""


def test_names_do_not_clash_with_structure_tables(tables):
    """Display names are free; the SQL names are prefixed."""
    table = tables.create("atom", columns=[("x", "float", 1.0)])
    assert table.sql_name.startswith("table_")
    table.append_row(x=2.0)
    assert tables.system_db["atom"] is not None
    assert tables["atom"].get_cell(table.current_row, "x") == 2.0


def test_create_existing(tables):
    tables.create("T")
    with pytest.raises(KeyError):
        tables.create("T")
    table = tables.create("T", columns=[("a", "integer", None)], replace=True)
    assert table.columns == ["a"]
    assert tables.names == ["T"]


def test_quoting(tables):
    """Column names with quotes, spaces and units work."""
    name = 'He said "hi" (kJ/mol)'
    table = tables.create("Odd 'name'", columns=[(name, "string", "it's")])
    rowid = table.append_row()
    assert table.get_cell(rowid, name) == "it's"


def test_append_defaults_and_current(tables):
    table = tables.create(
        "T",
        columns=[
            ("b", "boolean", None),
            ("i", "integer", None),
            ("f", "float", None),
            ("s", "string", None),
        ],
    )
    first, second = table.append_rows([{"i": 3}, {"s": "x", "b": True}])
    assert table.current_row == second
    assert table.get_row(first) == {
        "b": False,
        "i": 3,
        "f": table.get_cell(first, "f"),
        "s": "",
    }
    assert math.isnan(table.get_cell(first, "f"))
    assert table.get_cell(second, "b") is True
    assert table.get_cell(second, "i") == 0
    assert table.row_ids() == [first, second]
    assert table.position(second) == 1
    assert table.rowid_at(1) == second
    assert table.rowid_at(2) is None
    assert table.next_rowid(first) == second
    assert table.next_rowid(second) is None


def test_append_unknown_column(tables):
    table = tables.create("T", columns=[("a", "integer", None)])
    with pytest.raises(KeyError, match="no column"):
        table.append_row(b=1)


def test_add_column_fills_existing_rows(tables):
    table = tables.create("T", columns=[("a", "integer", None)])
    table.append_rows([{"a": 1}, {"a": 2}])
    assert table.add_column("c", "integer", 7)
    assert not table.add_column("c", "float", 1.0)
    assert [row["c"] for _, row in table.iter_rows()] == [7, 7]
    assert table.column_type("c") == "integer"


def test_set_get_cell(tables):
    table = tables.create("T", columns=[("a", "float", None), ("j", "json", None)])
    rowid = table.append_row()
    table.set_cell(rowid, "a", np.float64(1.5))
    assert table.get_cell(rowid, "a") == 1.5
    table.set_cell(rowid, "a", 2)
    assert table.get_cell(rowid, "a") == 2.0
    assert isinstance(table.get_cell(rowid, "a"), float)
    table.set_cell(rowid, "j", [1, 2, 3])
    assert table.get_cell(rowid, "j") == "[1,2,3]"
    with pytest.raises(KeyError):
        table.set_cell(rowid + 1, "a", 1.0)
    with pytest.raises(KeyError):
        table.set_cell(rowid, "nope", 1.0)


def test_text_is_kept_exactly(tables):
    """Untyped columns: text that looks like a number stays text."""
    table = tables.create("T", columns=[("v", "float", None)])
    rowid = table.append_row(v="1.0960")
    assert table.get_cell(rowid, "v") == "1.0960"
    df = table.to_dataframe()
    assert df["v"].dtype != "float64"
    assert df["v"].iloc[0] == "1.0960"


def test_numpy_values(tables):
    table = tables.create("T", columns=[("i", "integer", None), ("b", "boolean", None)])
    rowid = table.append_row(i=np.int64(4), b=np.bool_(True))
    assert table.get_cell(rowid, "i") == 4
    assert table.get_cell(rowid, "b") is True


def test_find_and_where(tables):
    table = tables.create("T", columns=[("k", "string", None), ("v", "integer", None)])
    ids = table.append_rows(
        [{"k": "a", "v": 1}, {"k": "b", "v": 2}, {"k": "a", "v": 3}]
    )
    assert table.find("k", "a") == [ids[0], ids[2]]
    assert [r["v"] for _, r in table.iter_rows(where=[("v", ">=", 2)])] == [2, 3]
    with pytest.raises(ValueError):
        list(table.iter_rows(where=[("v", "like", 2)]))


def test_index_column(tables):
    table = tables.create("T", columns=[("k", "string", None), ("v", "integer", None)])
    with pytest.raises(ValueError):
        table.index_column = "nope"
    table.index_column = "k"
    table.append_rows([{"k": "x", "v": 1}, {"k": "y", "v": 2}])
    df = table.to_dataframe()
    assert df.index.name == "k"
    assert list(df.index) == ["x", "y"]
    assert list(df.columns) == ["v"]
    assert df["v"].dtype == "int64"


def test_to_dataframe_types(tables):
    table = tables.create(
        "T",
        columns=[
            ("b", "boolean", None),
            ("i", "integer", None),
            ("f", "float", None),
            ("s", "string", None),
        ],
    )
    table.append_rows([{"b": True, "i": 1, "f": 1.5, "s": "x"}, {}])
    df = table.to_dataframe()
    assert list(df.columns) == ["b", "i", "f", "s"]
    assert df["b"].dtype == bool
    assert df["i"].dtype == "int64"
    assert df["f"].dtype == "float64"
    assert list(df.index) == [0, 1]
    assert df.to_csv(index=False) == "b,i,f,s\nTrue,1,1.5,x\nFalse,0,,\n"


def test_integer_with_null_becomes_float(tables):
    table = tables.create("T", columns=[("i", "integer", None)])
    rowid = table.append_row(i=1)
    table.append_row()
    table.set_cell(rowid, "i", None)
    assert table.to_dataframe()["i"].dtype == "float64"


def test_empty_dataframe(tables):
    table = tables.create("T")
    table.append_rows([{}, {}])
    df = table.to_dataframe()
    assert len(df) == 2
    assert list(df.columns) == []


def test_from_dataframe(tables):
    df = pandas.DataFrame(
        {"name": ["a", "b"], "n": [1, 2], "x": [0.5, math.nan], "ok": [True, False]}
    )
    table = tables.from_dataframe("T", df, index_column="name")
    assert [table.column_type(c) for c in table.columns] == [
        "string",
        "integer",
        "float",
        "boolean",
    ]
    assert table.current_row == table.rowid_at(0)
    out = table.to_dataframe()
    pandas.testing.assert_frame_equal(out, df.set_index("name"))


def test_current_row(tables):
    table = tables.create("T", columns=[("a", "integer", None)])
    rowid = table.append_row(a=1)
    table.current_row = None
    assert table.current_row is None
    table.current_row = rowid
    with pytest.raises(KeyError):
        table.current_row = rowid + 10


def test_metadata(tables):
    table = tables.create("T", metadata={"filename": "t.csv"})
    assert table.metadata == {"filename": "t.csv"}
    table.set_metadata("loop index", True)
    table.set_metadata("filename", None)
    assert table.metadata == {"loop index": True}


def test_journal(tables):
    table = tables.create("T", columns=[("a", "integer", None)])
    rowid = table.append_row(a=1)
    table.set_cell(rowid, "a", 2)
    table.add_column("b")
    ops = [(e["op"], e["row"], e["column"]) for e in tables.journal("T")]
    assert ops == [
        ("create", None, None),
        ("add_column", None, "a"),
        ("append", rowid, None),
        ("set", rowid, "a"),
        ("add_column", None, "b"),
    ]


def test_delete(tables):
    table = tables.create("T")
    sql_name = table.sql_name
    tables.delete("T")
    assert "T" not in tables
    assert sql_name not in tables.system_db
    assert tables.journal("T")[-1]["op"] == "drop"


def test_persistence(tmp_path):
    """Tables survive closing and reopening the file."""
    path = tmp_path / "seamm.db"
    db = SystemDB(filename=str(path))
    table = db.user_tables.create("T", columns=[("a", "integer", None)])
    table.append_row(a=5)
    db.db.commit()
    db.close()

    db = SystemDB(filename=str(path))
    table = db.user_tables["T"]
    assert table.to_dataframe()["a"].tolist() == [5]
    assert table.current_row == table.rowid_at(0)
    db.close()


def test_read_only(tmp_path):
    """A read-only database can be read but not written."""
    path = tmp_path / "seamm.db"
    db = SystemDB(filename=str(path))
    db.user_tables.create("T", columns=[("a", "integer", None)]).append_row(a=1)
    db.db.commit()
    db.close()

    db = SystemDB(filename=f"file:{path}?mode=ro")
    tables = db.user_tables
    assert tables.read_only
    assert tables["T"].to_dataframe()["a"].tolist() == [1]
    with pytest.raises(PermissionError, match="read-only"):
        tables["T"].append_row(a=2)
    with pytest.raises(PermissionError, match="read-only"):
        tables.create("U")
    db.close()


def test_columns_differing_in_case(tables):
    """Display names are case sensitive, as pandas columns were."""
    table = tables.create("T", columns=[("E", "float", None)])
    assert table.add_column("e", "float", None)
    rowid = table.append_row(E=1.0, e=2.0)
    assert table.get_cell(rowid, "E") == 1.0
    assert table.get_cell(rowid, "e") == 2.0
    assert list(table.to_dataframe().columns) == ["E", "e"]


def test_text_columns_store_text(tables):
    """An integer written to a text column is found by its text."""
    table = tables.create(
        "T",
        columns=[("name", "string", None), ("x", "float", None)],
        index_column="name",
    )
    rowid = table.append_row(name=3, x=1.0)
    table.append_row(name=True)
    assert table.find("name", "3") == [rowid]
    assert table.find("name", 3) == [rowid]
    assert table.to_dataframe().index.tolist() == ["3", "True"]


def test_missing_text_is_nan(tables):
    table = tables.create("T", columns=[("s", "string", None)])
    rowid = table.append_row()
    table.set_cell(rowid, "s", None)
    assert table.get_cell(rowid, "s") is None
    df = table.to_dataframe()
    assert math.isnan(df["s"].iloc[0])
    assert "NaN" in df.to_string()


def test_bad_index_column_changes_nothing(tables):
    """Validation happens before an existing table is replaced."""
    tables.create("T", columns=[("a", "integer", None)]).append_row(a=1)
    with pytest.raises(ValueError, match="index column"):
        tables.create(
            "T", columns=[("b", "integer", None)], index_column="c", replace=True
        )
    with pytest.raises(ValueError, match="Column type"):
        tables.create("U", columns=[("b", "decimal", None)])
    assert tables.names == ["T"]
    assert tables["T"].to_dataframe()["a"].tolist() == [1]


def test_row_ids_not_reused(tables):
    table = tables.create("T", columns=[("a", "integer", None)])
    first, second = table.append_rows([{"a": 1}, {"a": 2}])
    table.db.execute(f"DELETE FROM {table.table} WHERE __rowid__ = ?", (second,))
    assert table.append_row(a=3) == second + 1


def test_list_default(tables):
    table = tables.create("T", columns=[("j", "json", [1, 2])])
    rowid = table.append_row()
    assert table.get_cell(rowid, "j") == "[1,2]"


def test_from_dataframe_definitions(tables):
    df = pandas.DataFrame({"j": ["[1]", "[2]"], "n": [1, 2]})
    table = tables.from_dataframe(
        "T", df, definitions=[{"name": "j", "type": "json", "default": ""}]
    )
    assert table.column_type("j") == "json"
    assert table.column_type("n") == "integer"

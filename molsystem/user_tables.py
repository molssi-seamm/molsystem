# -*- coding: utf-8 -*-

"""User tables: named, ordered tables of results stored in the database.

These are the tables a flowchart builds (with the Table step, or with
``store_results``), as opposed to the tables that hold the structures. They are
kept in the same SQLite file, under SQL names ``table_<n>``, with a registry
``_tables`` holding each table's display name, declared column types and
defaults, index column, current row and free-form metadata, and a journal
``_table_changes`` recording every change.

The columns have no SQL type, so SQLite stores exactly what was written (the
text "1.0960" stays text); the declared type in the registry decides how values
are read back. Rows are identified by an internal row id, which is not meant to
be shown to users: they see the index column or the position of the row.

Writes are not committed here; the caller (the flowchart evaluator) commits
after each step.
"""

import json
import logging
import math
from typing import Any, Iterable, Iterator

import numpy as np
import pandas

from .table import _Table

logger = logging.getLogger(__name__)

#: The declared column types and the default for a column without one.
column_types = {
    "boolean": False,
    "integer": 0,
    "float": math.nan,
    "string": "",
    "json": "",
}

ROWID = "__rowid__"
REGISTRY = "_tables"
JOURNAL = "_table_changes"


def quote(name: str) -> str:
    """Quote an SQL identifier."""
    return '"' + str(name).replace('"', '""') + '"'


def to_sql(value: Any) -> Any:
    """Convert a Python value into one that sqlite3 can store."""
    if value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if isinstance(value, (bool, int, str, bytes)):
        return value
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, separators=(",", ":"))
    if isinstance(value, np.ndarray):
        return json.dumps(value.tolist(), separators=(",", ":"))
    try:
        if pandas.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def from_sql(value: Any, coltype: str) -> Any:
    """Convert a stored value back according to the declared column type."""
    if coltype == "boolean":
        if isinstance(value, int) and value in (0, 1):
            return bool(value)
    elif coltype == "float":
        if value is None:
            return math.nan
        if isinstance(value, int) and not isinstance(value, bool):
            return float(value)
    return value


def _cast_column(values: list, coltype: str) -> pandas.Series:
    """A pandas Series for a column, typed as its declared type allows."""
    if coltype == "float":
        if all(isinstance(v, (int, float)) for v in values):
            return pandas.Series(values, dtype="float64")
    elif coltype == "integer":
        if all(isinstance(v, int) and not isinstance(v, bool) for v in values):
            return pandas.Series(values, dtype="int64")
        if all(v is None or isinstance(v, (int, float)) for v in values):
            return pandas.Series(
                [math.nan if v is None else v for v in values], dtype="float64"
            )
    elif coltype == "boolean":
        if all(isinstance(v, bool) for v in values):
            return pandas.Series(values, dtype="bool")
    if len(values) == 0:
        return pandas.Series(values, dtype="object")
    # Let pandas infer, as it would have for the in-memory tables (e.g. str).
    return pandas.Series(values)


class UserTables:
    """The user tables in a SystemDB, by display name.

    Access through ``SystemDB.user_tables``. The registry and journal are created
    when the first table is.
    """

    def __init__(self, system_db):
        self._system_db = system_db
        self._items = {}

    def __contains__(self, name: str) -> bool:
        return self._has_registry() and self._exists(name)

    def __getitem__(self, name: str) -> "UserTable":
        if name not in self:
            raise KeyError(f"There is no table '{name}'.")
        if name not in self._items:
            self._items[name] = UserTable(self, name)
        return self._items[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self.names)

    def __len__(self) -> int:
        return len(self.names)

    @property
    def system_db(self):
        return self._system_db

    @property
    def db(self):
        return self._system_db.db

    @property
    def cursor(self):
        return self._system_db.cursor

    @property
    def names(self) -> list:
        """The names of the tables, in the order they were created."""
        if not self._has_registry():
            return []
        return [
            row[0]
            for row in self.db.execute(
                f"SELECT name FROM {REGISTRY} ORDER BY CAST(substr(sql_name, 7) AS"
                " INTEGER)"
            )
        ]

    @property
    def read_only(self) -> bool:
        """Whether the database is read-only, so tables cannot be written."""
        filename = self._system_db.filename
        return filename is not None and "mode=ro" in filename

    def check_writable(self, name: str = None):
        """Raise a clear error if tables cannot be written."""
        if self.read_only:
            what = "tables" if name is None else f"the table '{name}'"
            raise PermissionError(
                f"Cannot write {what}: the database {self._system_db.filename} was "
                "opened read-only, and tables are stored in the database."
            )

    def create(
        self,
        name: str,
        columns: Iterable = (),
        index_column: str = None,
        metadata: dict = None,
        replace: bool = False,
    ) -> "UserTable":
        """Create a new, empty table.

        Parameters
        ----------
        name : str
            The display name of the table.
        columns : [(name, type, default)] or [dict]
            The columns, in order. The default may be None for the type's default.
        index_column : str = None
            The column whose values identify the rows, if any.
        metadata : dict = None
            Free-form JSON-serializable information kept with the table.
        replace : bool = False
            Replace an existing table of the same name rather than raising.
        """
        self.check_writable(name)
        self._ensure_registry()
        if name in self:
            if not replace:
                raise KeyError(f"The table '{name}' already exists.")
            self.delete(name)

        n = self.db.execute(
            f"SELECT MAX(CAST(substr(sql_name, 7) AS INTEGER)) FROM {REGISTRY}"
        ).fetchone()[0]
        n = 1 if n is None else n + 1
        while f"table_{n}" in self._system_db:
            n += 1
        sql_name = f"table_{n}"

        self.db.execute(
            f"CREATE TABLE {quote(sql_name)} ({quote(ROWID)} INTEGER PRIMARY KEY)"
        )
        self.db.execute(
            f"INSERT INTO {REGISTRY} (name, sql_name, index_column, current_row,"
            " columns, metadata) VALUES (?, ?, ?, NULL, '[]', ?)",
            (name, sql_name, None, json.dumps({} if metadata is None else metadata)),
        )
        self._journal(name, "create")
        table = UserTable(self, name)
        self._items[name] = table
        for column in columns:
            if isinstance(column, dict):
                table.add_column(
                    column["name"], column.get("type", "string"), column.get("default")
                )
            else:
                table.add_column(*column)
        if index_column is not None:
            table.index_column = index_column
        return table

    def delete(self, name: str):
        """Remove a table."""
        self.check_writable(name)
        if name not in self:
            raise KeyError(f"There is no table '{name}'.")
        sql_name = self._sql_name(name)
        self.db.execute(f"DROP TABLE {quote(sql_name)}")
        self.db.execute(f"DELETE FROM {REGISTRY} WHERE name = ?", (name,))
        self._journal(name, "drop")
        self._items.pop(name, None)

    def from_dataframe(
        self,
        name: str,
        df: pandas.DataFrame,
        index_column: str = None,
        metadata: dict = None,
        replace: bool = False,
    ) -> "UserTable":
        """Create a table holding the contents of a DataFrame.

        A named index of the DataFrame becomes a column; ``index_column`` names the
        column that identifies rows.
        """
        if df.index.name is not None:
            df = df.reset_index()
        columns = []
        for column, dtype in zip(df.columns, df.dtypes):
            coltype = dtype_to_type(dtype)
            columns.append((str(column), coltype, None))
        table = self.create(
            name,
            columns=columns,
            index_column=index_column,
            metadata=metadata,
            replace=replace,
        )
        names = [str(c) for c in df.columns]
        rows = [dict(zip(names, values)) for values in df.itertuples(index=False)]
        table.append_rows(rows, move_current=False)
        table.current_row = table.rowid_at(0) if table.n_rows > 0 else None
        return table

    # Internal helpers
    def _has_registry(self) -> bool:
        return REGISTRY in self._system_db

    def _exists(self, name: str) -> bool:
        row = self.db.execute(
            f"SELECT COUNT(*) FROM {REGISTRY} WHERE name = ?", (name,)
        ).fetchone()
        return row[0] == 1

    def _ensure_registry(self):
        if self._has_registry():
            return
        self.db.execute(
            f"CREATE TABLE {REGISTRY} ("
            "  name TEXT PRIMARY KEY,"
            "  sql_name TEXT UNIQUE NOT NULL,"
            "  index_column TEXT,"
            "  current_row INTEGER,"
            "  columns TEXT NOT NULL DEFAULT '[]',"
            "  metadata TEXT NOT NULL DEFAULT '{}'"
            ")"
        )
        self.db.execute(
            f"CREATE TABLE IF NOT EXISTS {JOURNAL} ("
            "  seq INTEGER PRIMARY KEY AUTOINCREMENT,"
            '  "table" TEXT NOT NULL,'
            '  "row" INTEGER,'
            '  "column" TEXT,'
            '  "op" TEXT NOT NULL'
            ")"
        )

    def _sql_name(self, name: str) -> str:
        row = self.db.execute(
            f"SELECT sql_name FROM {REGISTRY} WHERE name = ?", (name,)
        ).fetchone()
        if row is None:
            raise KeyError(f"There is no table '{name}'.")
        return row[0]

    def _journal(self, name, op, row=None, column=None):
        self.db.execute(
            f'INSERT INTO {JOURNAL} ("table", "row", "column", "op")'
            " VALUES (?, ?, ?, ?)",
            (name, row, column, op),
        )

    def journal(self, name: str = None) -> list:
        """The journal entries as dicts, optionally for one table, in order."""
        if not self._has_registry():
            return []
        sql = f'SELECT seq, "table", "row", "column", "op" FROM {JOURNAL}'
        parameters = ()
        if name is not None:
            sql += ' WHERE "table" = ?'
            parameters = (name,)
        sql += " ORDER BY seq"
        return [
            dict(zip(("seq", "table", "row", "column", "op"), row))
            for row in self.db.execute(sql, parameters)
        ]


def dtype_to_type(dtype) -> str:
    """The declared column type for a pandas dtype."""
    kind = getattr(dtype, "kind", "O")
    if kind == "b":
        return "boolean"
    if kind in "iu":
        return "integer"
    if kind == "f":
        return "float"
    return "string"


class UserTable(_Table):
    """One user table. Get it from ``SystemDB.user_tables[name]``."""

    def __init__(self, tables: UserTables, name: str):
        self._tables = tables
        self._name = name
        super().__init__(tables.system_db, tables._sql_name(name))

    def __repr__(self) -> str:
        return f"UserTable('{self._name}', {self.n_rows} rows)"

    __str__ = __repr__

    @property
    def name(self) -> str:
        """The display name of the table."""
        return self._name

    @property
    def sql_name(self) -> str:
        return self._table

    # Registry-backed properties
    def _registry(self, field):
        return self.db.execute(
            f"SELECT {field} FROM {REGISTRY} WHERE name = ?", (self._name,)
        ).fetchone()[0]

    def _set_registry(self, field, value):
        self._tables.check_writable(self._name)
        self.db.execute(
            f"UPDATE {REGISTRY} SET {field} = ? WHERE name = ?", (value, self._name)
        )

    @property
    def column_definitions(self) -> list:
        """The columns as an ordered list of dicts: name, type, default."""
        return json.loads(self._registry("columns"))

    @property
    def columns(self) -> list:
        """The column names, in order."""
        return [c["name"] for c in self.column_definitions]

    @property
    def defaults(self) -> dict:
        return {c["name"]: c["default"] for c in self.column_definitions}

    def column_type(self, column: str) -> str:
        for c in self.column_definitions:
            if c["name"] == column:
                return c["type"]
        raise KeyError(f"The table '{self._name}' has no column '{column}'.")

    @property
    def index_column(self):
        return self._registry("index_column")

    @index_column.setter
    def index_column(self, value):
        if value is not None and value not in self.columns:
            columns = ", ".join(self.columns)
            raise ValueError(
                f"The index column '{value}' is not in the table '{self._name}': "
                f"columns = {columns}"
            )
        self._set_registry("index_column", value)

    @property
    def current_row(self):
        """The row id of the current row, or None for 'the next row to append'."""
        return self._registry("current_row")

    @current_row.setter
    def current_row(self, rowid):
        if rowid is not None and not self.has_row(rowid):
            raise KeyError(f"The table '{self._name}' has no row with id {rowid}.")
        self._set_registry("current_row", rowid)

    @property
    def metadata(self) -> dict:
        """A copy of the table's metadata. Change it with set_metadata."""
        return json.loads(self._registry("metadata"))

    def set_metadata(self, key: str, value):
        metadata = self.metadata
        if value is None:
            metadata.pop(key, None)
        else:
            metadata[key] = value
        self._set_registry("metadata", json.dumps(metadata))

    # Columns
    def add_column(self, name: str, coltype: str = "string", default=None) -> bool:
        """Add a column, filling existing rows with its default.

        Returns False, changing nothing, if the column already exists.
        """
        name = str(name)
        if coltype not in column_types:
            raise ValueError(
                f"Column type '{coltype}' must be one of {', '.join(column_types)}."
            )
        definitions = self.column_definitions
        if any(c["name"] == name for c in definitions):
            return False
        if name == ROWID:
            raise ValueError(f"'{ROWID}' is reserved and cannot name a column.")
        self._tables.check_writable(self._name)
        if default is None:
            default = column_types[coltype]
        if isinstance(default, float) and math.isnan(default):
            stored_default = None
        else:
            stored_default = default
        self.db.execute(f"ALTER TABLE {self.table} ADD COLUMN {quote(name)}")
        if stored_default is not None:
            self.db.execute(
                f"UPDATE {self.table} SET {quote(name)} = ?", (to_sql(stored_default),)
            )
        definitions.append({"name": name, "type": coltype, "default": stored_default})
        self._set_registry("columns", json.dumps(definitions))
        self._tables._journal(self._name, "add_column", column=name)
        return True

    def default(self, column: str):
        """The default value of a column, as it would be read back."""
        for c in self.column_definitions:
            if c["name"] == column:
                return from_sql(c["default"], c["type"])
        raise KeyError(f"The table '{self._name}' has no column '{column}'.")

    # Rows
    def row_ids(self) -> list:
        """The row ids, in order."""
        return [
            row[0]
            for row in self.db.execute(
                f"SELECT {quote(ROWID)} FROM {self.table} ORDER BY {quote(ROWID)}"
            )
        ]

    def has_row(self, rowid) -> bool:
        return (
            self.db.execute(
                f"SELECT COUNT(*) FROM {self.table} WHERE {quote(ROWID)} = ?",
                (rowid,),
            ).fetchone()[0]
            == 1
        )

    def rowid_at(self, position: int):
        """The row id at a 0-based position, or None past the end."""
        if position < 0:
            return None
        row = self.db.execute(
            f"SELECT {quote(ROWID)} FROM {self.table} ORDER BY {quote(ROWID)}"
            " LIMIT 1 OFFSET ?",
            (position,),
        ).fetchone()
        return None if row is None else row[0]

    def position(self, rowid) -> int:
        """The 0-based position of a row."""
        if not self.has_row(rowid):
            raise KeyError(f"The table '{self._name}' has no row with id {rowid}.")
        return self.db.execute(
            f"SELECT COUNT(*) FROM {self.table} WHERE {quote(ROWID)} < ?", (rowid,)
        ).fetchone()[0]

    def find(self, column: str, value) -> list:
        """The row ids whose column equals the value, in order."""
        self.column_type(column)
        return [
            row[0]
            for row in self.db.execute(
                f"SELECT {quote(ROWID)} FROM {self.table} WHERE {quote(column)} = ?"
                f" ORDER BY {quote(ROWID)}",
                (to_sql(value),),
            )
        ]

    def next_rowid(self, rowid):
        """The row id after the given one, or None at the end."""
        row = self.db.execute(
            f"SELECT MIN({quote(ROWID)}) FROM {self.table} WHERE {quote(ROWID)} > ?",
            (rowid,),
        ).fetchone()
        return None if row is None else row[0]

    def append_rows(self, rows: Iterable[dict], move_current: bool = True) -> list:
        """Append rows given as dicts; missing columns get their defaults.

        Returns the new row ids. By default the current row becomes the last one.
        """
        self._tables.check_writable(self._name)
        definitions = self.column_definitions
        names = [c["name"] for c in definitions]
        defaults = [c["default"] for c in definitions]
        known = set(names)
        parameters = []
        for row in rows:
            unknown = [k for k in row if k not in known]
            if len(unknown) > 0:
                raise KeyError(
                    f"The table '{self._name}' has no column(s) "
                    f"{', '.join(repr(k) for k in unknown)}. The columns are: "
                    f"{', '.join(names)}"
                )
            parameters.append(
                [
                    to_sql(row[name]) if name in row else default
                    for name, default in zip(names, defaults)
                ]
            )
        if len(parameters) == 0:
            return []
        last = self.db.execute(f"SELECT MAX({quote(ROWID)}) FROM {self.table}")
        last = last.fetchone()[0]
        first = 1 if last is None else last + 1
        rowids = [*range(first, first + len(parameters))]
        columns = ", ".join(quote(c) for c in [ROWID, *names])
        places = ", ".join(["?"] * (len(names) + 1))
        self.db.executemany(
            f"INSERT INTO {self.table} ({columns}) VALUES ({places})",
            [[rowid, *p] for rowid, p in zip(rowids, parameters)],
        )
        self.db.executemany(
            f'INSERT INTO {JOURNAL} ("table", "row", "column", "op")'
            " VALUES (?, ?, NULL, 'append')",
            [(self._name, rowid) for rowid in rowids],
        )
        if move_current:
            self._set_registry("current_row", rowids[-1])
        return rowids

    def append_row(self, move_current: bool = True, **values):
        """Append one row; returns its row id."""
        return self.append_rows([values], move_current=move_current)[0]

    def set_cell(self, rowid, column: str, value):
        """Set one value."""
        self._tables.check_writable(self._name)
        self.column_type(column)
        cursor = self.db.execute(
            f"UPDATE {self.table} SET {quote(column)} = ? WHERE {quote(ROWID)} = ?",
            (to_sql(value), rowid),
        )
        if cursor.rowcount != 1:
            raise KeyError(f"The table '{self._name}' has no row with id {rowid}.")
        self._tables._journal(self._name, "set", row=rowid, column=column)

    def get_cell(self, rowid, column: str):
        """Get one value, converted according to the column's declared type."""
        coltype = self.column_type(column)
        row = self.db.execute(
            f"SELECT {quote(column)} FROM {self.table} WHERE {quote(ROWID)} = ?",
            (rowid,),
        ).fetchone()
        if row is None:
            raise KeyError(f"The table '{self._name}' has no row with id {rowid}.")
        return from_sql(row[0], coltype)

    def get_row(self, rowid) -> dict:
        for _rowid, row in self.iter_rows(where=[(ROWID, "=", rowid)]):
            return row
        raise KeyError(f"The table '{self._name}' has no row with id {rowid}.")

    def iter_rows(self, where: Iterable = ()) -> Iterator:
        """Iterate over (row id, dict of values) in order.

        ``where`` is a list of (column, op, value) with op one of the SQL
        comparisons =, !=, <, <=, >, >=, all of which must hold.
        """
        definitions = self.column_definitions
        names = [c["name"] for c in definitions]
        types = [c["type"] for c in definitions]
        sql = f"SELECT {', '.join(quote(c) for c in [ROWID, *names])} FROM {self.table}"
        clauses = []
        parameters = []
        for column, op, value in where:
            if op == "==":
                op = "="
            if op not in ("=", "!=", "<", "<=", ">", ">="):
                raise ValueError(f"Unsupported comparison '{op}'.")
            clauses.append(f"{quote(column)} {op} ?")
            parameters.append(to_sql(value))
        if len(clauses) > 0:
            sql += " WHERE " + " AND ".join(clauses)
        sql += f" ORDER BY {quote(ROWID)}"
        for line in self.db.execute(sql, parameters).fetchall():
            yield line[0], {
                name: from_sql(value, coltype)
                for name, coltype, value in zip(names, types, line[1:])
            }

    def to_dataframe(self) -> pandas.DataFrame:
        """A copy of the table as a DataFrame, typed by the declared types.

        The index is the index column if there is one, otherwise 0, 1, 2, ...
        """
        definitions = self.column_definitions
        data = {c["name"]: [] for c in definitions}
        for _rowid, row in self.iter_rows():
            for name, value in row.items():
                data[name].append(value)
        df = pandas.DataFrame(
            {c["name"]: _cast_column(data[c["name"]], c["type"]) for c in definitions}
        )
        if len(definitions) == 0:
            df = pandas.DataFrame(index=pandas.RangeIndex(self.n_rows))
        index = self.index_column
        if index is not None:
            df.set_index(index, inplace=True)
        return df

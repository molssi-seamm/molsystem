# -*- coding: utf-8 -*-

import collections.abc
from itertools import zip_longest
import logging
import pandas

from .column import _Column

logger = logging.getLogger(__name__)

column_types = {"int": "INTEGER", "float": "REAL", "str": "TEXT", "bytes": "BLOB"}


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Table(collections.abc.MutableMapping):
    """A dictionary-like object for holding tabular data

    This is a wrapper around an SQL table, with columns being the attributes.

    Attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be deleted. See the method 'add_attribute' for more detail.

    Rows can be added ('append') or deleted ('delete').
    """

    def __init__(self, system_db, table: str, other=None, logger=logger):
        self._system_db = system_db
        self.logger = logger
        if "." in table:
            self._schema, self._table = table.split(".")
        else:
            self._schema = "main"
            self._table = table
        if other is not None:
            self.copy(other)
        self._checkpoints = []

    def __enter__(self):
        """Copy the table to a backup for a 'with' statement."""
        n = len(self._checkpoints)
        checkpoint = f"temp.{self._table}_checkpoint_{n}"
        self.db.execute(f"CREATE TEMP TABLE {checkpoint} AS SELECT * FROM {self.table}")
        self.db.commit()
        self._checkpoints.append(checkpoint)
        return self

    def __exit__(self, etype, value, traceback) -> None:
        """Handle returning from a 'with' statement."""
        checkpoint = self._checkpoints.pop()

        if etype is None:
            self.db.commit()

            # Log the changes
            previous = _Table(self.system_db, checkpoint)
            diffs = self.diff(previous)
            if diffs is not None:
                self.logger.debug(diffs)
        else:
            # Delete everything in the table and copy the checkpoint back.
            self.db.execute("PRAGMA foreign_keys = OFF")
            self.clear()
            self.db.execute(f"INSERT INTO {self.table} SELECT * FROM {checkpoint}")
            self.db.execute("PRAGMA foreign_keys = ON")
            self.db.commit()

        # Drop the temporary table
        self.db.execute(f"DROP TABLE {checkpoint}")
        self.db.commit()

        return False

    def __getitem__(self, key):
        """Allow [] to access the data!"""
        return self.get_column(key)

    def __setitem__(self, key, value) -> None:
        """Allow x[key] access to the data"""
        self[key][0:] = value

    def __delitem__(self, key) -> None:
        """Allow deletion of keys"""
        # The easy way, which is not supported in SQLite :-(
        # self.cursor.execute(f'ALTER TABLE {self.table} DROP {key}')

        columns = set(self.attributes)
        columns.remove(key)
        column_def = ", ".join(columns)

        # Need the unquoted version!
        table = self._table

        sql = f"""
        -- disable foreign key constraint check
        PRAGMA foreign_keys=off;

        -- Here you can drop column
        CREATE TABLE {self.schema}."tmp_{table}"
           AS SELECT {column_def} FROM {self.table};

        -- drop the table
        DROP TABLE {self.table};

        -- rename the new_table to the table
        ALTER TABLE {self.schema}."tmp_{table}"
          RENAME TO "{self._table}";

        -- enable foreign key constraint check
        PRAGMA foreign_keys=on;
        """
        self.db.executescript(sql)
        self.db.commit()

    def __iter__(self) -> iter:
        """Allow iteration over the object"""
        return iter([*self.attributes.keys()])

    def __len__(self) -> int:
        """The len() command"""
        return len(self.attributes)

    def __repr__(self) -> str:
        """The string representation of this object"""
        df = self.to_dataframe()
        return repr(df)

    def __str__(self) -> str:
        """The pretty string representation of this object"""
        df = self.to_dataframe()
        return str(df)

    def __contains__(self, item) -> bool:
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.

        tmp_item = item.strip('"')
        return tmp_item in self.attributes

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        # Check numbers of rows
        if self.n_rows != other.n_rows:
            return False

        # Check the columns
        attributes = self.attributes
        columns = set(attributes)

        other_attributes = other.attributes
        other_columns = set(other_attributes)

        if columns != other_columns:
            return False

        is_same = True
        # Need to check the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        db = self.system_db
        other_db = other.system_db

        detach = False
        if db.filename != other_db.filename:
            if db.is_attached(other_db):
                schema = db.attached_as(other_db)
            else:
                # Attach the other system_db in order to do comparisons.
                schema = self.system_db.attach(other_db)
                detach = True
            other_table = f'"{schema}"."{other._table}"'
        else:
            other_table = other.table
        table = self.table

        self.cursor.execute(
            f"""
            SELECT COUNT(*) FROM
            (SELECT rowid, * FROM {other_table}
            EXCEPT
            SELECT rowid, * FROM {table})
            """
        )
        if self.cursor.fetchone()[0] != 0:
            is_same = False
        else:
            self.cursor.execute(
                f"""
                SELECT COUNT(*) FROM
                (SELECT rowid, * FROM {table}
                EXCEPT
                SELECT rowid, * FROM {other_table})
                """
            )
            if self.cursor.fetchone()[0] != 0:
                is_same = False

        # Detach the other database if needed
        if detach:
            self.system_db.detach(other_db)

        return is_same

    @property
    def system_db(self):
        """The system_db that we belong to."""
        return self._system_db

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def cursor(self):
        """The database connection."""
        return self.system_db.cursor

    @property
    def schema(self):
        """The database schema the table is in."""
        return f'"{self._schema}"'

    @property
    def table(self) -> str:
        """The name of this table."""
        return f'"{self._schema}"."{self._table}"'

    @property
    def n_rows(self) -> int:
        """The number of rows in the table."""
        self.cursor.execute(f"SELECT COUNT(*) FROM {self.table}")
        result = self.cursor.fetchone()[0]
        return result

    @property
    def attributes(self):
        """The definitions of the attributes."""
        result = {}
        for row in self.db.execute(
            "   SELECT *"
            f"    FROM {self.schema}.pragma_table_info('{self._table}')"
            " ORDER BY 'cid'"
        ):
            result[row["name"]] = {
                "type": row["type"],
                "notnull": bool(row["notnull"]),
                "default": row["dflt_value"],
                "primary key": bool(row["pk"]),
            }
        for row in self.db.execute(f"PRAGMA foreign_key_list('{self.table}')"):
            if row["to"] is None:
                result[row["from"]]["fk"] = row["table"]
            else:
                result[row["from"]]["fk"] = f"{row['table']}.{row['to']}"

        return result

    # @property
    # def version(self):
    #     return self.system_db.version

    def add_attribute(
        self,
        name: str,
        coltype: str = "float",
        default=None,
        notnull: bool = False,
        index: bool = False,
        pk: bool = False,
        references: str = None,
        on_delete: str = "cascade",
        on_update: str = "cascade",
        values=None,
    ) -> None:
        """Adds a new attribute.

        If the default value is None, you must always provide values wherever
        needed, for example when adding a row.

        Parameters
        ----------
        name : str
            The name of the attribute.
        coltype : str
            The type of the attribute (column). Must be one of 'int',
            'float', 'str' or 'byte'
        default : int, float, str or byte
            The default value for the attribute if no value is given.
        notnull : bool = False
            Whether the value must be non-null
        index : bool = False
            Whether to create an index on the column
        pk : bool = False
            Whether the column is the primry keys
        references : str = None
            If not null, the column is a foreign key for this table.
        on_delete : str = 'cascade'
            How to handle deletions of a foregin keys
        on_update : str = 'cascade'
            How to handle updates of a foregin key
        values : Any
            Either a single value or a list of values length 'nrows' of
            values to fill the column.

        Returns:
            None
        """
        # Does the table exist?
        if self.table in self.system_db:
            table_exists = True
            if pk:
                raise ValueError(
                    "The primary key can only be set on the first attribute "
                    "of a table that does not yet exist."
                )
        else:
            table_exists = False
            if pk:
                notnull = True

        if name in self:
            raise RuntimeError(
                "{} attribute '{}' is already defined!".format(
                    self.__class__.__name__, name
                )
            )

        # not null columns must have defaults
        if not pk and notnull and default is None:
            raise ValueError(f"Not null attributes must have defaults: {name}")

        # see if the values are given
        if values is not None:
            length = self.length_of_values(values)

            if length > 1 and length != self.n_rows:
                raise IndexError(
                    f"The number of values given, {length}, must be either 1, "
                    f"or the number of rows in {self.table}: {self.n_rows}"
                )

        # Create the column
        parameters = []
        if coltype in column_types:
            column_type = column_types[coltype]
        else:
            column_type = coltype
        column_def = f'"{name}" {column_type}'
        if default is not None:
            if coltype == "str":
                column_def += f" DEFAULT '{default}'"
            else:
                column_def += f" DEFAULT {default}"
        if notnull:
            column_def += " NOT NULL"
        if references is not None:
            column_def += f" REFERENCES {references}"
            if on_delete is not None and on_delete != "":
                column_def += f" ON DELETE {on_delete}"
            if on_update is not None and on_update != "":
                column_def += f" ON UPDATE {on_update}"

        if table_exists:
            self.cursor.execute(
                f"ALTER TABLE {self.table} ADD {column_def}", parameters
            )
        else:
            if pk:
                column_def += " PRIMARY KEY"
            self.cursor.execute(f"CREATE TABLE {self.table} ({column_def})", parameters)

        if index == "unique":
            self.cursor.execute(
                f'CREATE UNIQUE INDEX {self.schema}."idx_{name}"'
                f'  ON "{self._table}"'
                f"  ('{name}')"
            )
        elif index:
            self.cursor.execute(
                f'CREATE INDEX {self.schema}."idx_{name}"'
                f'  ON "{self._table}"'
                f"  ('{name}')"
            )

        if values is not None:
            self[name] = values

        self.db.commit()

    def append(self, n=None, **kwargs):
        """Append one or more rows

        The keywords are the names of attributes and the value to use.
        The default value for any attributes not given is used unless it is
        'None' in which case an error is thrown. It is an error if there is not
        an exisiting attribute corresponding to any given as arguments.

        Parameters
        ----------
        kwargs : keyword arguments
            Giving existing attributes and values.

        Returns
        -------
        [int]
            The ids of the created rows.
        """

        if len(kwargs) > 0:
            n_rows, lengths = self._get_n_rows(**kwargs)
        else:
            lengths = {}
            if n is not None:
                n_rows = n
            else:
                n_rows = 1

        if n is not None:
            if n_rows != 1 and n_rows != n:
                raise RuntimeError(
                    f"Requested number of rows ({n}) not compatible with the "
                    f"length of the data ({n_rows})."
                )
            n_rows = n

        # Check that any missing attributes have defaults
        attributes = self.attributes
        for key in attributes:
            if not attributes[key]["notnull"] or attributes[key]["primary key"]:
                continue

            if key not in kwargs:
                if (
                    "default" not in attributes[key]
                    or attributes[key]["default"] is None
                ):
                    raise KeyError(
                        "There is no default for attribute "
                        "'{}'. You must supply a value".format(key)
                    )

        # Add id's if needed
        if "id" not in kwargs and "id" in attributes:
            self.cursor.execute(f"SELECT MAX(id) FROM {self.table}")
            last_id = self.cursor.fetchone()[0]
            if last_id is None:  # Table is empty
                last_id = 0
            kwargs["id"] = [*range(last_id + 1, last_id + n_rows + 1)]
            lengths["id"] = n_rows

        # All okay, so proceed.
        values = {}
        for key, value in kwargs.items():
            if lengths[key] == 0:
                values[key] = [value] * n_rows
            elif lengths[key] == 1:
                values[key] = [value[0]] * n_rows
            else:
                values[key] = value

        parameters = []
        for row in range(n_rows):
            line = []
            for key, value in values.items():
                line.append(value[row])
            parameters.append(line)

        columns = '"' + '", "'.join(kwargs.keys()) + '"'
        places = ", ".join(["?"] * len(values.keys()))

        self.cursor.executemany(
            f"INSERT INTO {self.table} ({columns}) VALUES ({places})", parameters
        )

        self.db.commit()

        if "id" in kwargs:
            return kwargs["id"]

    def clear(self):
        """Delete all rows from the table."""
        self.db.execute(f"DELETE FROM {self.table}")

    def diff(self, other):
        """Difference between this table and another

        Parameters
        ----------
        other : _Table
            The other table to diff against

        Result
        ------
        result : dict
            The differences, described in a dictionary
        """
        self.logger.debug("diffs for {self.table}")
        result = {}
        summary = {}

        # Check the columns
        attributes = self.attributes
        columns = set(attributes)

        other_attributes = other.attributes
        other_columns = set(other_attributes)

        if columns == other_columns:
            column_def = "rowid, *"
        else:
            added = columns - other_columns
            if len(added) > 0:
                self.logger.debug(f"      columns added: {', '.join(sorted(added))}")
                result["columns added"] = added
                summary["columns add"] = added
            deleted = other_columns - columns
            if len(deleted) > 0:
                self.logger.debug(f"     columns deleted: {', '.join(sorted(deleted))}")
                result["columns deleted"] = deleted
                summary["columns deleted"] = deleted

            in_common = other_columns & columns
            if len(in_common) > 0:
                column_def = "rowid, " + ", ".join([f'"{x}"' for x in in_common])
            else:
                # No columns shared
                return result

        # Need to check the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        db = self.system_db
        other_db = other.system_db

        detach = False
        if db.filename != other_db.filename:
            if db.is_attached(other_db):
                schema = db.attached_as(other_db)
            else:
                # Attach the other system_db in order to do comparisons.
                schema = self.system_db.attach(other_db)
                detach = True
            other_table = f'"{schema}"."{other._table}"'
        else:
            other_table = other.table
        table = self.table

        changed = {}
        last = None
        sql = f"""
        SELECT {column_def} FROM
         (
          SELECT {column_def} FROM {other_table}
           EXCEPT
          SELECT {column_def} FROM {table}
         )
        UNION ALL
        SELECT {column_def} FROM
         (
          SELECT {column_def} FROM {table}
           EXCEPT
          SELECT {column_def} FROM {other_table}
         )
        ORDER BY rowid
        """
        changed_columns = set()
        for row in self.db.execute(sql):
            if last is None:
                last = row
            elif row["rowid"] == last["rowid"]:
                changes = set()
                for k1, v1, v2 in zip(last.keys(), last, row):
                    if v1 != v2:
                        changed_columns.add(k1)
                        changes.add((k1, v1, v2))
                changed[row["rowid"]] = changes
                last = None
            else:
                last = row
        if len(changed) > 0:
            self.logger.debug(
                f"    columns changed: {', '.join(sorted(changed_columns))}"
            )
            result["changed"] = changed
            summary["columns changed"] = changed_columns

        # See about the rows added
        added = {}
        if "id" in self:
            for row in self.db.execute(
                f"""
                SELECT * FROM {table}
                WHERE rowid NOT IN (SELECT rowid FROM {other_table})
                """
            ):
                added[row["id"]] = row[1:]
        else:
            for row in self.db.execute(
                f"""
                SELECT rowid, * FROM {table}
                WHERE rowid NOT IN (SELECT rowid FROM {other_table})
                """
            ):
                added[row["rowid"]] = row[1:]

        if len(added) > 0:
            self.logger.debug(f"         rows added: {len(added)}")
            result["columns in added rows"] = row.keys()[1:]
            result["added"] = added
            summary["rows added"] = len(added)

        # See about the rows deleted
        deleted = {}
        if "id" in self:
            for row in self.db.execute(
                f"""
                SELECT * FROM {other_table}
                WHERE rowid NOT IN (SELECT rowid FROM {table})
                """
            ):
                deleted[row["id"]] = row[1:]
        else:
            for row in self.db.execute(
                f"""
                SELECT rowid, * FROM {other_table}
                WHERE rowid NOT IN (SELECT rowid FROM {table})
                """
            ):
                deleted[row["rowid"]] = row[1:]

        if len(deleted) > 0:
            self.logger.debug(f"       rows deleted: {len(deleted)}")
            result["columns in deleted rows"] = row.keys()[1:]
            result["deleted"] = deleted
            summary["rows deleted"] = len(deleted)

        # Detach the other database if needed
        if detach:
            self.system_db.detach(other_db)

        if len(result) == 0:
            self.logger.debug("    no changes")
            result = None
        else:
            result["summary"] = summary

        return result

    def get_as_dict(self, *args):
        """Return the table data as a Python dictionary of lists.

        Parameters
        ----------
        args : [str]
            Added selection criteria for the SQL, one word at a time.

        Returns
        -------
        dict(str: [])
            A dictionary whose keys are the column names and values as lists
        """
        rows = self.rows(*args)
        columns = [x[0] for x in rows.description]
        data = {key: [] for key in columns}
        for row in rows:
            for key, value in zip(columns, row):
                data[key].append(value)

        return data

    def get_column(self, column, sql=None, where=""):
        """Return a column of data from the table.

        Parameters
        ----------
        column : str
            The name of the column.
        sql : str = None
            Optional SQL statement to get the rows.
        where : str = ''
            Optional WHERE clause when getting the rows if no SQL is given.
            Should be the complete clause, starting with 'WHERE ...'
        """
        return _Column(self, column=column, sql=sql, where=where)

    def get_column_data(self, column, sql=None, where=""):
        """Return a column of data from the table as a list.

        Parameters
        ----------
        column : str
            The name of the column.
        sql : str = None
            Optional SQL statement to get the rows.
        where : str = ''
            Optional WHERE clause when getting the rows if no SQL is given.
            Should be the complete clause, starting with 'WHERE ...'
        """
        result = []
        if sql is None:
            sql = f"SELECT {column} FROM {self.table} {where}"

        for row in self.db.execute(sql):
            result.append(row[0])
        return result

    def delete(self, *args):
        """Delete rows matching the selection."""
        if len(args) == 0:
            return self.db.execute(f"DELETE FROM {self.table}")

        sql = f"DELETE FROM {self.table} WHERE"

        parameters = []
        clause = []
        for col, op, value in grouped(args, 3):
            if op == "==":
                op = "="
            clause.append(f'"{col}" {op} ?')
            parameters.append(value)

        sql += " AND ".join(clause)
        return self.db.execute(sql, parameters)

    def rows(self, *args):
        """Return an iterator over the rows."""
        if len(args) == 0:
            return self.db.execute(f"SELECT * FROM {self.table}")

        sql = f"SELECT * FROM {self.table} WHERE"

        parameters = []
        clause = []
        for col, op, value in grouped(args, 3):
            if op == "==":
                op = "="
            clause.append(f'"{col}" {op} ?')
            parameters.append(value)

        sql += " AND ".join(clause)
        return self.db.execute(sql, parameters)

    def _get_n_rows(self, **kwargs):
        """Get the total number of rows represented in the arguments."""
        n_rows = None

        lengths = {}
        for key, value in kwargs.items():
            if key not in self:
                raise KeyError(
                    f'"{key}" is not an attribute of the table ' f"{self.table}!"
                )
            length = self.length_of_values(value)
            lengths[key] = length

            if n_rows is None:
                n_rows = 1 if length == 0 else length

            if length > 1 and length != n_rows:
                if n_rows == 1:
                    n_rows = length
                else:
                    raise IndexError(
                        'key "{}" has the wrong number of values, '.format(key)
                        + "{}. Should be 1 or the number of rows in {} ({}).".format(
                            length, self.table, n_rows
                        )
                    )
        return n_rows, lengths

    def copy(self, other):
        """Replace this object with a copy of another."""
        if self.table in self.system_db:
            del self.system_db[self.table]

        # Need the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.

        name = self.system_db.nickname
        other_name = other.system_db.nickname
        detach = False
        if name != other_name:
            if not self.system_db.is_attached(other_name):
                # Attach the other system_db in order to do comparisons.
                self.system_db.attach(other.system_db)
                detach = True
            other_table = f'"{other_name}"."{other._table}"'
        else:
            other_table = other.table
        table = self.table

        self.cursor.execute(f"CREATE TABLE {table} AS SELECT * FROM {other_table}")
        self.db.commit()

        # Detach the other database if needed
        if detach:
            self.system_db.detach(other.system_db)

    def length_of_values(self, values):
        """Return the length of the values argument.

        Parameters
        ----------
        values : Any
            The values, which might be a string, single value, list, tuple,
            etc.

        Returns
        -------
        length : int
            The length of the values. 0 indicates a scalar
        """
        if isinstance(values, str):
            return 0
        else:
            try:
                return len(values)
            except TypeError:
                return 0

    def to_dataframe(self):
        """Return the contents of the table as a Pandas Dataframe."""
        data = {}
        for line in self.cursor.execute(f"SELECT rowid, * FROM {self.table}"):
            data[line[0]] = line[1:]
        columns = [x[0] for x in self.cursor.description[1:]]

        df = pandas.DataFrame.from_dict(data, orient="index", columns=columns)

        return df

# -*- coding: utf-8 -*-

import csv
import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


class _Properties(object):
    """A class for handling the properties in the star schema."""

    def __init__(self, system_db):
        self._system_db = system_db
        self._standard_properties = None

    @property
    def cursor(self):
        return self.system_db.cursor

    @property
    def db(self):
        return self.system_db.db

    @property
    def standard_properties(self):
        """The standard properties recognized by SEAMM.

        These are officially defined properties that can be used anywhere in SEAMM, as
        long as the type and definition correspond to the standard.

        You can create other properties on te fly, but they must be prefixed by a unique
        name followed by a dot ('.') so that they do not conflict with either the
        standard properties or other properties defined on-the-fly. Typically the unique
        name should  that of the program generating the property, which implicitly
        defines details of the property.

        For example, the standard property "enthalpy of formation" refers to the
        experimental heat of formation, or a calculated value comparable to experimental
        values. If you are not sure what the heat of formation in e.g. MOPAC is, you
        could create a new property "MOPAC.enthalpy of formation", which is clearly
        similar to the standard "enthalpy of formation". If the community decides that
        it is indeed the same, it can be replaced by the standard form, and also aliased
        to it for backwards compatibility.
        """
        if self._standard_properties is None:
            self._standard_properties = {}
            path = Path(pkg_resources.resource_filename(__name__, "data/"))
            csv_file = path / "standard_properties.csv"
            with open(csv_file, newline="", encoding="utf-8-sig") as fd:
                data = csv.reader(fd)
                line = 0
                for row in data:
                    line += 1
                    if line == 1:
                        # Headers
                        headers = [*row]
                        if headers != [
                            "Property",
                            "Type",
                            "Units",
                            "Description",
                            "URL",
                        ]:
                            raise ValueError(
                                "Header of standard properties file not valid: "
                                + ", ".join(headers)
                            )
                    else:
                        property = row[0]
                        data = self._standard_properties[property] = {}
                        for key, value in zip(headers[1:], row[1:]):
                            data[key] = value
        return self._standard_properties

    @property
    def system_db(self):
        """Return the SystemDB object."""
        return self._system_db

    def add(self, name, _type="float", units=None, description="", noerror=False):
        """Add a property to the database. By default, it is an error if it already
        exists!

        Parameters
        ----------
        name : str
            The name of the property.
        _type : str = "float"
            The type of the property -- float, int, str. Default is "float".
        units : str = None
            The units of the property value, as a string that Pint can interpret.
        description : str
            A longer, text description of the property.
        noerror: bool = False
            Whether to quietly ignore an existing entry

        Returns
        -------
        int
            The database id of the property.

        Note
        ----
        If the property is registered in the list of standard properties, just
        give the name and the rest will be filled in from the standard property
        metadata.
        """
        if self.exists(name):
            if noerror:
                return self.property_id(name)
            else:
                raise ValueError(f"Property '{name}' already exists.")

        table = self.system_db["property"]
        if name in self.standard_properties:
            data = self.standard_properties[name]
            result = table.append(
                name=name,
                type=data["Type"],
                units=data["Units"],
                description=data["Description"],
            )
        else:
            result = table.append(
                name=name,
                type=_type,
                units=units,
                description=description,
            )
        return result[0]

    def create_schema(self):
        """Add the needed tables to the database."""
        # The property dimension table
        table = self.system_db["property"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("name", coltype="str", index="unique")
        table.add_attribute("type", coltype="str")
        table.add_attribute("units", coltype="str")
        table.add_attribute("description", coltype="str")

        # Floating point facts
        table = self.system_db["float_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="float")

        self.db.execute(
            "CREATE INDEX float_data_idx_configuration_property_value"
            "    ON float_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_property_value"
            "    ON float_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_configuration_property"
            "    ON float_data(configuration, property)"
        )

        # Integer facts
        table = self.system_db["int_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="int")

        self.db.execute(
            "CREATE INDEX int_data_idx_configuration_property_value"
            "    ON int_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_property_value"
            "    ON int_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_configuration_property"
            "    ON int_data(configuration, property)"
        )

        # String facts
        table = self.system_db["str_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="str")

        self.db.execute(
            "CREATE INDEX str_data_idx_configuration_property_value"
            "    ON str_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_property_value"
            "    ON str_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_configuration_property"
            "    ON str_data(configuration, property)"
        )

    def exists(self, name):
        """Whether the named property exists.

        Parameters
        ----------
        name : str
            The name of the property.

        Returns
        -------
        bool
            True if the property is registered, False otherwise.
        """
        self.cursor.execute("SELECT COUNT(*) FROM property WHERE name = ?", (name,))
        return self.cursor.fetchone()[0] != 0

    def get(self, configuration_id, _property):
        """Get the given property value for the configuration.

        Parameters
        ----------
        configuration_id : int
            The id of the configuration.
        _property : int or str
            The id or name of the property.

        Returns
        -------
        int, float, or str
            The value of the property.
        """
        if isinstance(_property, str):
            if self.exists(_property):
                pid = self.property_id(_property)
            else:
                raise ValueError(f"Property '{_property}' does not exist.")
        else:
            pid = _property

        ptype = self.property_type(pid)
        sql = (
            f"SELECT value FROM {ptype}_data"
            "  WHERE configuration = ? AND property = ?"
        )
        self.cursor.execute(sql, (configuration_id, pid))
        result = self.cursor.fetchone()
        if result is None:
            raise ValueError(
                f"Property {_property} does not exist for configuration "
                f"{configuration_id}"
            )
        return result[0]

    def known_properties(self):
        """List the known properties."""
        self.cursor.execute("SELECT name FROM property")
        return [row[0] for row in self.cursor.fetchall()]

    def property_id(self, name):
        """The id for a property

        Parameters
        ----------
        name : str
            The name of the property.

        Returns
        -------
        int
            The database id for the property.
        """
        self.cursor.execute("SELECT id FROM property WHERE name = ?", (name,))
        return self.cursor.fetchone()[0]

    def property_name(self, pid):
        """The name of a property

        Parameters
        ----------
        pid : int
            The id of the property.

        Returns
        -------
        str
            The name of the property.
        """
        self.cursor.execute("SELECT name FROM property WHERE id = ?", (pid,))
        return self.cursor.fetchone()[0]

    def property_type(self, _property):
        """The type of a property

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        str
            The type of the property.
        """
        if isinstance(_property, str):
            self.cursor.execute(
                "SELECT type FROM property WHERE name = ?", (_property,)
            )
        else:
            self.cursor.execute("SELECT type FROM property WHERE id = ?", (_property,))
        return self.cursor.fetchone()[0]

    def put(self, configuration_id, _property, value):
        """Store the given property value for the configuration.

        Parameters
        ----------
        configuration_id : int
            The id of the configuration.
        _property : int or str
            The id or name of the property.
        value : int, float, or str
            The value to store.
        """
        if isinstance(_property, str):
            if self.exists(_property):
                pid = self.property_id(_property)
            else:
                raise ValueError(f"Property '{_property}' does not exist.")
        else:
            pid = _property

        ptype = self.property_type(pid)
        sql = (
            f"INSERT INTO {ptype}_data (configuration, property, value)"
            "      VALUES(?, ?, ?)"
        )
        self.db.execute(sql, (configuration_id, pid, value))

    def query(self, *args, what=["configuration"]):
        """Find configurations that match the query defined by the args.

        what : [str or int]]
            What to return ... "configuration", property names or ids.
        where : {}

        Note
        ----
        args should be a multiple of 3 in length, with each triplet being
        (property, operator, value). The property
        """
        sql = "SELECT"
        where = "WHERE"
        criteria = []
        tables = {}
        results = []
        n_tables = 0
        for i, item in enumerate(what):
            if i > 0:
                sql += ","
            if isinstance(item, str):
                if item == "configuration":
                    sql += " t0.configuration"
                    results.append("configuration_id")
                else:
                    pid = self.property_id(item)
                    ptype = self.property_type(pid)
                    table = ptype + "_data"
                    if table not in tables:
                        tables[table] = {pid: f"t{n_tables}"}
                        n_tables += 1
                    elif pid not in tables[table]:
                        tables[table][pid] = f"t{n_tables}"
                        n_tables += 1
                    alias = tables[table][pid]
                    sql += f" {alias}.value"
                    criteria.append(f"{alias}.property == {pid}")
                    results.append(item)
            elif isinstance(item, int):
                pid = item
                ptype = self.property_type(pid)
                table = ptype + "_data"
                if table not in tables:
                    tables[table] = {pid: f"t{n_tables}"}
                    n_tables += 1
                elif pid not in tables[table]:
                    tables[table][pid] = f"t{n_tables}"
                    n_tables += 1
                alias = tables[table][pid]
                sql += f" {alias}.value"
                criteria.append(f"{alias}.property == {pid}")
                results.append(self.property_name(pid))

        # And the WHERE clause ...
        where = {0: ""}
        items = iter(args)
        level = 0
        for item in items:
            if item == "(":
                level += 1
                where[level] = " ("
            elif item == ")":
                where[level] += ")"
                where[level - 1] += where[level]
                del where[level]
                level -= 1
            else:
                # Configuration or Property name or id
                if isinstance(item, str):
                    pid = self.property_id(item)
                    ptype = self.property_type(pid)
                    table = ptype + "_data"
                    if table not in tables:
                        tables[table] = {pid: f"t{n_tables}"}
                        n_tables += 1
                    elif pid not in tables[table]:
                        tables[table][pid] = f"t{n_tables}"
                        n_tables += 1
                    alias = tables[table][pid]
                    where[level] += f" {alias}.property == {pid}"
                elif isinstance(item, int):
                    pid = item
                    ptype = self.property_type(pid)
                    table = ptype + "_data"
                    if table not in tables:
                        tables[table] = {pid: f"t{n_tables}"}
                        n_tables += 1
                    elif pid not in tables[table]:
                        tables[table][pid] = f"t{n_tables}"
                        n_tables += 1
                    alias = tables[table][pid]
                    where[level] += f" {alias}.property == {pid}"

                operator = next(items)
                if operator.lower() == "between":
                    if ptype == "str":
                        where[level] += (
                            f' AND {alias}.value BETWEEN "{next(items)}"'
                            f' AND "{next(items)}"'
                        )
                    else:
                        where[level] += (
                            f" AND {alias}.value BETWEEN {next(items)}"
                            f" AND {next(items)}"
                        )
                else:
                    if ptype == "str":
                        where[level] += f' AND {alias}.value {operator} "{next(items)}"'
                    else:
                        where[level] += f" AND {alias}.value {operator} {next(items)}"

        # Put it all together ... FROM clause plus needed criteria
        # if len(tables) == 0:
        #     # No tables ... so probably getting all configurations with properties
        #     tables = ["int_data", "float_data", "str_data"]
        sql += " FROM "

        tmp = []
        for table, data in tables.items():
            for pid, alias in data.items():
                tmp.append(f"{table} AS {alias}")
        sql += ", ".join(tmp)

        # Now for the WHERE clause
        sql += " WHERE "
        tmp = []
        for table, data in tables.items():
            for pid, alias in data.items():
                if alias != "t0":
                    tmp.append(f"{alias}.configuration = t0.configuration")
        sql += " AND ".join(tmp)

        if len(tmp) == 0:
            sql += " "
        else:
            sql += " AND "
        sql += " AND ".join(criteria)

        # The user criteria, if any
        if where[0] != "":
            sql += " AND"
            sql += where[0]

        print(sql)

        result = {item: [] for item in results}
        self.cursor.execute(sql)
        for row in self.cursor:
            for item, value in zip(results, row):
                result[item].append(value)
        return result

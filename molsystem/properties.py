# -*- coding: utf-8 -*-

import csv
import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401

logger = logging.getLogger(__name__)
standard_properties = {}


def add_properties_from_file(path):
    """The standard properties recognized by SEAMM.

    These are officially defined properties that can be used anywhere in SEAMM, as
    long as the type and definition correspond to the standard.

    Each property is defined by a string with up to three parts:

        <property name>#<code or 'experiment'>#<technique or model chemistry>

    The property name is required. In most cases this is followed by either 'experiment'
    or the name of the code, e.g. 'MOPAC', 'Gaussian', or 'VASP'. The final part,
    if present, is either the experimental technique used to measure the property, or
    the model chemistry, such as 'MP2/6-31G**', 'PM7', or a forcefield name such as
    'AMBER/ff19SB'.

    You can create other properties on the fly, but they follow the above convention
    and should have an appropriate code and, if necessary, model chemistry, so that
    they full name is unique and does not conflict with any other defined name.

    For example, the standard property "enthalpy of formation" refers to the
    experimental heat of formation, or a calculated value comparable to experimental
    values. If you are not sure what the heat of formation in e.g. MOPAC is, you could
    create a new property "enthalpy of formation#MOPAC#<parameterization>", which is
    clearly similar to the standard "enthalpy of formation". If the community decides
    that it is indeed the same, it can be replaced by the standard form, and also
    aliased to it for backwards compatibility.
    """
    global standard_properties
    with open(path, newline="", encoding="utf-8-sig") as fd:
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
                data = standard_properties[property] = {}
                for key, value in zip(headers[1:], row[1:]):
                    data[key] = value


path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "standard_properties.csv"
add_properties_from_file(csv_file)


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
        global standard_properties
        return standard_properties

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
        table.add_attribute("system", coltype="int", references="system")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="float")

        self.db.execute(
            "CREATE INDEX float_data_idx_configuration_property_value"
            "    ON float_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_system_property_value"
            "    ON float_data(system, property, value)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_property_value"
            "    ON float_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_configuration_property"
            "    ON float_data(configuration, property)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_system_property"
            "    ON float_data(system, property)"
        )
        self.db.execute(
            "CREATE INDEX float_data_idx_configuration ON float_data(configuration)"
        )
        self.db.execute("CREATE INDEX float_data_idx_system ON float_data(system)")

        # Integer facts
        table = self.system_db["int_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("system", coltype="int", references="system")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="int")

        self.db.execute(
            "CREATE INDEX int_data_idx_configuration_property_value"
            "    ON int_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_system_property_value"
            "    ON int_data(system, property, value)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_property_value"
            "    ON int_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_configuration_property"
            "    ON int_data(configuration, property)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_system_property ON int_data(system, property)"
        )
        self.db.execute(
            "CREATE INDEX int_data_idx_configuration ON int_data(configuration)"
        )
        self.db.execute("CREATE INDEX int_data_idx_system ON int_data(system)")

        # String facts
        table = self.system_db["str_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("system", coltype="int", references="system")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="str")

        self.db.execute(
            "CREATE INDEX str_data_idx_configuration_property_value"
            "    ON str_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_system_property_value"
            "    ON str_data(system, property, value)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_property_value"
            "    ON str_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_configuration_property"
            "    ON str_data(configuration, property)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_system_property  ON str_data(system, property)"
        )
        self.db.execute(
            "CREATE INDEX str_data_idx_configuration ON str_data(configuration)"
        )
        self.db.execute("CREATE INDEX str_data_idx_system ON str_data(system)")

    def description(self, _property):
        """The description of a property

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        str
            The description of the property.
        """
        return self.metadata(_property)[2]

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

    def get(self, configuration_id, _property="all"):
        """Get the given property value for the configuration.

        Parameters
        ----------
        configuration_id : int
            The id of the configuration.
        _property : int or str, or "all"
            The id or name of the property, or all properties if "all".

        Returns
        -------
        int, float, or str
            The value of the property.
        """
        if _property == "all":
            sql = "SELECT name, type, value"
            sql += "  FROM property, float_data"
            sql += " WHERE float_data.property = property.id AND configuration = ?"
            sql += " UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, int_data"
            sql += " WHERE int_data.property = property.id AND configuration = ?"
            sql += " UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, str_data"
            sql += " WHERE str_data.property = property.id AND configuration = ?"

            self.cursor.execute(
                sql, (configuration_id, configuration_id, configuration_id)
            )

            result = {}
            for row in self.cursor:
                name, _type, value = row
                if _type == "float":
                    result[name] = float(value)
                elif _type == "int":
                    result[name] = int(value)
                else:
                    result[name] = value
            return result
        else:
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

    def get_for_system(self, system_id, _property="all"):
        """Get the given property value(s) for the system.

        Parameters
        ----------
        system_id : int
            The id of the system.
        _property : int or str of "all"
            The id or name of the property, or "all" for all properties.

        Returns
        -------
        [int, float, or str]
            The value(s) of the property.
        """
        if _property == "all":
            sql = "SELECT name, type, value"
            sql += "  FROM property, float_data"
            sql += " WHERE float_data.property = property.id AND system = ?"
            sql += " UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, int_data"
            sql += " WHERE int_data.property = property.id AND system = ?"
            sql += " UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, str_data"
            sql += " WHERE str_data.property = property.id AND system = ?"

            self.cursor.execute(sql, (system_id, system_id, system_id))

            result = {}
            for row in self.cursor:
                name, _type, value = row
                if _type == "float":
                    result[name] = float(value)
                elif _type == "int":
                    result[name] = int(value)
                else:
                    result[name] = value
        else:
            if isinstance(_property, str):
                if self.exists(_property):
                    pid = self.property_id(_property)
                else:
                    raise ValueError(f"Property '{_property}' does not exist.")
            else:
                pid = _property
            ptype = self.property_type(pid)

            sql = f"SELECT value FROM {ptype}_data WHERE system = ? AND property = ?"
            result = []
            for row in self.db.execute(sql, (system_id, pid)):
                result.append(row[0])

        return result

    def id(self, name):
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
        tmp = self.cursor.fetchone()
        if tmp is None:
            raise KeyError(f"Property '{name}' is not known.")
        else:
            return tmp[0]

    def known_properties(self):
        """List the known properties."""
        self.cursor.execute("SELECT name FROM property")
        return [row[0] for row in self.cursor.fetchall()]

    def metadata(self, _property):
        """The metadata for a property

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        str, str, str
            The type, units, and description of the property
        """
        if isinstance(_property, str):
            self.cursor.execute(
                "SELECT type, units, description FROM property WHERE name = ?",
                (_property,),
            )
        else:
            self.cursor.execute(
                "SELECT type, units, description FROM property WHERE id = ?",
                (_property,),
            )
        tmp = self.cursor.fetchone()
        if tmp is not None:
            return tmp

        if _property in self.standard_properties:
            data = self.standard_properties[_property]
            return data["Type"], data["Units"], data["Description"]

        raise KeyError(f"Property '{_property}' is not known.")

    def name(self, pid):
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
        tmp = self.cursor.fetchone()
        if tmp is None:
            raise KeyError(f"Property id = '{pid}' is not known.")
        else:
            return tmp[0]

    def property_id(self, name):
        "Obsolete routine kept for compatibility"
        return self.id(name)

    def property_name(self, pid):
        "Obsolete routine kept for compatibility"
        return self.name(pid)

    def property_type(self, _property):
        "Obsolete routine kept for compatibility"
        return self.type(_property)

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
        # Get the property id and type
        if isinstance(_property, str):
            if not self.exists(_property):
                if _property in self.standard_properties:
                    self.add(_property)
                else:
                    raise ValueError(f"Property '{_property}' does not exist.")
            pid = self.property_id(_property)
        else:
            pid = _property
        ptype = self.property_type(pid)

        # Get the system id
        self.cursor.execute(
            "SELECT system FROM configuration WHERE id = ?", (configuration_id,)
        )
        system_id = self.cursor.fetchone()[0]

        sql = (
            f"INSERT INTO {ptype}_data (configuration, system, property, value)"
            "      VALUES(?, ?, ?, ?)"
        )
        self.db.execute(sql, (configuration_id, system_id, pid, value))

    def put_for_system(self, system_id, _property, value):
        """Store the given property value for the system.

        Parameters
        ----------
        system_id : int
            The id of the system.
        _property : int or str
            The id or name of the property.
        value : int, float, or str
            The value to store.
        """
        # Get the property id and type
        if isinstance(_property, str):
            if not self.exists(_property):
                if _property in self.standard_properties:
                    self.add(_property)
                else:
                    raise ValueError(f"Property '{_property}' does not exist.")
            pid = self.property_id(_property)
        else:
            pid = _property
        ptype = self.property_type(pid)

        sql = (
            f"INSERT INTO {ptype}_data (system, property, value)"
            "      VALUES(?, ?, ?)"
        )
        self.db.execute(sql, (system_id, pid, value))

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
                    if where[level] != "" and where[level] != "( ":
                        where[level] += " AND"
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
                    if where[level] != "" and where[level] != "( ":
                        where[level] += " AND"
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

        if len(criteria) > 0:
            if len(tmp) > 0:
                sql += " AND "
            sql += " AND ".join(criteria)

        # The user criteria, if any
        if where[0] != "":
            if len(tmp) > 0:
                sql += " AND"
            sql += where[0]

        logger.info(sql)

        result = {item: [] for item in results}
        try:
            self.cursor.execute(sql)
        except Exception:
            print("")
            print(f"{tmp=}")
            print(f"{criteria=}")
            print(f"{where[0]=}")
            raise

        for row in self.cursor:
            for item, value in zip(results, row):
                result[item].append(value)
        return result

    def type(self, _property):
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
        return self.metadata(_property)[0]

    def units(self, _property):
        """The unit string of a property

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        str
            The units string for the property.
        """
        return self.metadata(_property)[1]

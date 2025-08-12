# -*- coding: utf-8 -*-

import csv
import json
import logging
from pathlib import Path
import pkg_resources

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
                return self.id(name)
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

        # Array facts
        table = self.system_db["json_data"]
        table.add_attribute("id", coltype="int", pk=True)
        table.add_attribute("configuration", coltype="int", references="configuration")
        table.add_attribute("system", coltype="int", references="system")
        table.add_attribute("property", coltype="int", references="property")
        table.add_attribute("value", coltype="str")

        self.db.execute(
            "CREATE INDEX json_data_idx_configuration_property_value"
            "    ON json_data(configuration, property, value)"
        )
        self.db.execute(
            "CREATE INDEX json_data_idx_system_property_value"
            "    ON json_data(system, property, value)"
        )
        self.db.execute(
            "CREATE INDEX json_data_idx_property_value"
            "    ON json_data(property, value)"
        )
        self.db.execute(
            "CREATE INDEX json_data_idx_configuration_property"
            "    ON json_data(configuration, property)"
        )
        self.db.execute(
            "CREATE INDEX json_data_idx_system_property ON json_data(system, property)"
        )
        self.db.execute(
            "CREATE INDEX json_data_idx_configuration ON json_data(configuration)"
        )
        self.db.execute("CREATE INDEX json_data_idx_system ON json_data(system)")

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

    def get(
        self,
        _id,
        pattern="*",
        match="glob",
        is_system=False,
        include_system_properties=False,
        include_configuration_properties=False,
        types=["float", "int", "str", "json"],
    ):
        """Get the property value(s)

        Parameters
        ----------
        _id : int
            The id of the configuration or system.
        is_system : bool=False
            Whether the ID is for a system, not a configuration.
        pattern : str = "*"
            The pattern of the property.
        match : str="glob"
            Whether to use exact, glob, or 'like' matching.
        include_system_properties : bool=False
            For a configuration, whether to include properties that are on the system.
        include_configuration_properties : bool=False
            For a system, whether to include properties that are on the configurations
            of the system.
        types : [str] = ["float", "int", "str", "json"]
            The type of results to return.

        Returns
        -------
        {str: {str: value}}
            The matching property values.
        """
        if match == "glob":
            op = "GLOB"
        elif match == "like":
            op = "LIKE"
        elif match == "exact":
            op = "="
        else:
            raise ValueError(f"Unknown match type '{match}'.")

        sql = ""
        args = []

        if is_system:
            # Handle properties of systems
            if include_configuration_properties:
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                    sql += (
                        "       SELECT name, type, value, system, configuration\n"
                        f"         FROM property, {_type}_data\n"
                        "        WHERE system = ?\n"
                        f"          AND {_type}_data.property = property.id\n"
                        f"          AND name {op} ?\n"
                    )
                    args.append(_id)
                    args.append(pattern)
            else:
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                    sql += (
                        "       SELECT name, type, value, system, configuration\n"
                        f"         FROM property, {_type}_data\n"
                        "        WHERE system = ?\n"
                        "          AND configuration is Null\n"
                        f"          AND {_type}_data.property = property.id\n"
                        f"          AND name {op} ?\n"
                    )
                    args.append(_id)
                    args.append(pattern)
        else:
            for _type in types:
                if len(args) > 1:
                    sql += "        UNION\n"
                sql += (
                    "       SELECT name, type, value, system, configuration\n"
                    f"         FROM property, {_type}_data\n"
                    "        WHERE configuration = ?\n"
                    f"          AND {_type}_data.property = property.id\n"
                    f"          AND name {op} ?\n"
                )
                args.append(_id)
                args.append(pattern)
            if include_system_properties:
                self.cursor.execute(
                    "SELECT system FROM configuration WHERE id = ?", (_id,)
                )
                system_id = self.cursor.fetchone()[0]
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                    sql += (
                        "       SELECT name, type, value, system, configuration\n"
                        f"         FROM property, {_type}_data\n"
                        "        WHERE system = ?\n"
                        "          AND configuration is Null\n"
                        f"          AND {_type}_data.property = property.id\n"
                        f"          AND name {op} ?\n"
                    )
                    args.append(system_id)
                    args.append(pattern)

        self.cursor.execute(sql, args)

        result = {}
        for row in self.cursor:
            name, _type, value, sid, cid = row
            result[name] = {
                "sid": sid,
                "cid": cid,
            }
            if _type == "float":
                try:
                    result[name]["value"] = float(value)
                except Exception as e:
                    logger.warning(f"Error with value of property '{name}': {str(e)}")
                    del result[name]
            elif _type == "int":
                try:
                    result[name]["value"] = int(value)
                except Exception as e:
                    logger.warning(f"Error with value of property '{name}': {str(e)}")
                    del result[name]
            elif _type == "json":
                try:
                    result[name]["value"] = json.loads(value)
                except Exception as e:
                    logger.warning(f"Error with value of property '{name}': {str(e)}")
                    del result[name]
            else:
                result[name]["value"] = value
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

    def list(
        self,
        _id,
        pattern="*",
        match="glob",
        is_system=False,
        include_system_properties=False,
        include_configuration_properties=False,
        as_ids=False,
        types=["float", "int", "str", "json"],
    ):
        """Get the list of matching properties for the configuration.

        Parameters
        ----------
        _id : int
            The id of the configuration.
        is_system : bool=False
            Whether the ID is for a system, not a configuration.
        pattern : str = "*"
            The pattern of the property.
        match : str="glob"
            Whether to use exact, glob, or 'like' matching.
        include_system_properties : bool=False
            For a configuration, whether to include properties that are on the system.
        include_configuration_properties : bool=False
            For a system, whether to include properties that are on the configurations
            of the system.
        as_ids : bool=False
            Whether to return the ids rather than names
        types : [str] = ["float", "int", "str", "json"]
            The type of results to return.

        Returns
        -------
        [str] or [int]
            The matching properties.
        """
        if match == "glob":
            op = "GLOB"
        elif match == "like":
            op = "LIKE"
        elif match == "exact":
            op = "="
        else:
            raise ValueError(f"Unknown match type '{match}'.")

        if as_ids:
            sql = "SELECT id FROM property\n"
        else:
            sql = "SELECT name FROM property\n"
        sql += f" WHERE name {op} ?\n" "    AND id IN (\n"

        args = [pattern]

        if is_system:
            # Handle properties of systems
            if include_configuration_properties:
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                        f"       SELECT property FROM {_type}_data WHERE system = ?\n"
                    args.append(_id)
            else:
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                    sql += (
                        "       SELECT property\n"
                        f"         FROM {_type}_data\n"
                        "        WHERE system = ? AND configuration is Null\n"
                    )
                    args.append(_id)
        else:
            for _type in types:
                if len(args) > 1:
                    sql += "        UNION\n"
                sql += f"  SELECT property FROM {_type}_data WHERE configuration = ?\n"
                args.append(_id)
            if include_system_properties:
                self.cursor.execute(
                    "SELECT system FROM configuration WHERE id = ?", (_id,)
                )
                system_id = self.cursor.fetchone()[0]
                for _type in types:
                    if len(args) > 1:
                        sql += "        UNION\n"
                    sql += (
                        "       SELECT property\n"
                        f"         FROM {_type}_data\n"
                        "        WHERE configuration is Null AND system = ?\n"
                    )
                    args.append(system_id)

        sql += "    )"

        self.cursor.execute(sql, args)

        return [row[0] for row in self.cursor]

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

    def put(self, _id, _property, value, is_system=False):
        """Store the given property value for the configuration or system

        If the property value already exists, overwrite it in situ.

        Parameters
        ----------
        _id : int
            The id of the configuration.
        _property : int or str
            The id or name of the property.
        value : int, float, or str
            The value to store.
        is_system : bool=False
            Whether the ID is for a system, not a configuration.
        """
        # Get the property id and type
        if isinstance(_property, str):
            if not self.exists(_property):
                if _property in self.standard_properties:
                    self.add(_property)
                else:
                    raise ValueError(f"Property '{_property}' does not exist.")
            pid = self.id(_property)
        else:
            pid = _property
        ptype = self.type(pid)

        if ptype == "json":
            value = json.dumps(value, separators=(",", ":"))

        if is_system:
            # Check if a value already exists
            sql = (
                f"SELECT COUNT(*) FROM {ptype}_data WHERE system = ? AND "
                "configuration is Null AND property = ?"
            )
            self.cursor.execute(sql, (_id, pid))
            if self.cursor.fetchone()[0] == 0:
                # Not there so add
                sql = (
                    f"INSERT INTO {ptype}_data (system, property, value)"
                    "      VALUES(?, ?, ?)"
                )
                self.db.execute(sql, (_id, pid, value))
            else:
                # update in place
                sql = (
                    f"UPDATE {ptype}_data SET value = ? WHERE system = ? AND "
                    "configuration is Null AND property = ?"
                )
                self.db.execute(sql, (value, _id, pid))
        else:
            # Get the system id
            self.cursor.execute("SELECT system FROM configuration WHERE id = ?", (_id,))
            sid = self.cursor.fetchone()[0]

            # Check if a value already exists
            sql = (
                f"SELECT COUNT(*) FROM {ptype}_data WHERE system = ? AND "
                "configuration = ? AND property = ?"
            )
            self.cursor.execute(sql, (sid, _id, pid))
            if self.cursor.fetchone()[0] == 0:
                # Not there so add
                sql = (
                    f"INSERT INTO {ptype}_data (configuration, system, property, value)"
                    "      VALUES(?, ?, ?, ?)"
                )
                self.db.execute(sql, (_id, sid, pid, value))
            else:
                # update in place
                sql = (
                    f"UPDATE {ptype}_data SET value = ? WHERE system = ? AND "
                    "configuration = ? AND property = ?"
                )
                self.db.execute(sql, (value, sid, _id, pid))

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
        types = []
        n_tables = 0
        for i, item in enumerate(what):
            if i > 0:
                sql += ","
            if isinstance(item, str):
                if item == "configuration":
                    sql += " t0.configuration"
                    results.append("configuration_id")
                    types.append["configuration"]
                else:
                    pid = self.id(item)
                    ptype = self.type(pid)
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
                    types.append(ptype)
            elif isinstance(item, int):
                pid = item
                ptype = self.type(pid)
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
                results.append(self.name(pid))
                types.append(ptype)

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
                    pid = self.id(item)
                else:
                    pid = item
                ptype = self.type(pid)
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
                    if ptype == "json":
                        value1 = json.dumps(next(items), separator=(",", ":"))
                        value2 = json.dumps(next(items), separator=(",", ":"))
                        where[
                            level
                        ] += f' AND {alias}.value BETWEEN "{value1}" AND "{value2}"'
                    else:
                        where[level] += (
                            f" AND {alias}.value BETWEEN {next(items)}"
                            f" AND {next(items)}"
                        )
                else:
                    if ptype == "str":
                        where[level] += f' AND {alias}.value {operator} "{next(items)}"'
                    elif ptype == "json":
                        value = json.dumps(next(items), separator=(",", ":"))
                        where[level] += f' AND {alias}.value {operator} "{value}"'
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
            for item, ptype, value in zip(results, types, row):
                if ptype == "json":
                    result[item].append(json.loads(value))
                else:
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

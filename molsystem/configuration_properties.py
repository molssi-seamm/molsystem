# -*- coding: utf-8 -*-

"""Property methods for configurations."""

import logging
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


class _ConfigurationProperties(object):
    """A class for handling the properties of a configuration in the star schema."""

    def __init__(self, configuration):
        self._configuration = configuration
        self._cid = configuration.id
        self._properties = None

    @property
    def configuration(self):
        """The configuration these properties are connected with."""
        return self._configuration

    @property
    def properties(self):
        """The general properties object."""
        if self._properties is None:
            self._properties = self.configuration.system_db.properties
        return self._properties

    @property
    def standard_properties(self):
        """A thin wrapper of the _Properties method."""
        return self.properties.standard_properties()

    def add(self, name, _type="float", units=None, description="", noerror=False):
        """A thin wrapper of the _Properties method."""
        return self.properties.add(
            name, _type=_type, units=units, description=description, noerror=noerror
        )

    def exists(self, name):
        """A thin wrapper of the _Properties method."""
        return self.properties.exists(name)

    def get(self, _property="all"):
        """Get the given property value for this configuration.

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        int, float, or str
            The value of the property.
        """
        if _property == "all":
            sql = "SELECT name, type, value"
            sql += "  FROM property, float_data"
            sql += " WHERE float_data.property = property.id and property.id IN ("
            sql += "     SELECT property FROM float_data WHERE configuration = ?"
            sql += ") UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, int_data"
            sql += " WHERE int_data.property = property.id and property.id IN ("
            sql += "     SELECT property FROM int_data WHERE configuration = ?"
            sql += ") UNION "
            sql += "SELECT name, type, value"
            sql += "  FROM property, str_data"
            sql += " WHERE str_data.property = property.id and property.id IN ("
            sql += "     SELECT property FROM str_data WHERE configuration = ?"
            sql += ")"

            self.properties.cursor.execute(sql, (self._cid, self._cid, self._cid))

            result = {}
            for row in self.properties.cursor:
                name, _type, value = row
                if _type == "float":
                    result[name] = float(value)
                elif _type == "int":
                    result[name] = int(value)
                else:
                    result[name] = value
            return result
        else:
            return self.properties.get(self._cid, _property)

    def known_properties(self):
        """List the known properties."""
        return self.properties.known_properties()

    def list(self):
        """Return all of the property names for the current configuration.
        Returns
        -------
        [str]
            The names of the properties for the configuration.
        """
        sql = "SELECT name FROM property WHERE id IN ("
        sql += "SELECT property FROM float_data WHERE configuration = ?"
        sql += "UNION SELECT property FROM int_data WHERE configuration = ?"
        sql += "UNION SELECT property FROM str_data WHERE configuration = ?"
        sql += ")"

        self.properties.cursor.execute(sql, (self._cid, self._cid, self._cid))
        return [row[0] for row in self.properties.cursor]

    def list_ids(self):
        """Return all of the property ids for the current configuration.
        Returns
        -------
        [int]
            The names of the properties for the configuration.
        """
        sql = "SELECT property FROM float_data WHERE configuration = ?"
        sql += "UNION SELECT property FROM int_data WHERE configuration = ?"
        sql += "UNION SELECT property FROM str_data WHERE configuration = ?"

        self.properties.cursor.execute(sql, (self._cid, self._cid, self._cid))
        return [row[0] for row in self.properties.cursor]

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
        return self.properties.property_id(name)

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
        return self.properties.property_name(pid)

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
        return self.properties.property_type(_property)

    def put(self, _property, value):
        """Store the given property value for this configuration.

        Parameters
        ----------
        _property : int or str
            The id or name of the property.
        value : int, float, or str
            The value to store.
        """
        self.properties.put(self._cid, _property, value)

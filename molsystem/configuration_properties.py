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
        return self.properties.description(_property)

    def exists(self, name):
        """A thin wrapper of the _Properties method."""
        return self.properties.exists(name)

    def get(self, _property="all", include_system_properties=False):
        """Get the given property value for this configuration.

        Parameters
        ----------
        _property : int or str
            The id or name of the property.
        include_system_properties : bool=False
            Whether to include properties that are on the system, not any configuration

        Returns
        -------
        int, float, or str
            The value of the property.
        """
        return self.properties.get(
            self._cid,
            _property,
            include_system_properties=include_system_properties,
        )

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
        return self.properties.id(name)

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

    def metadata(self, _property):
        """The metadata for a property

        Parameters
        ----------
        _property : int or str
            The id or name of the property.

        Returns
        -------
        str, str, str, str
            The name, type, units, and description of the property
        """
        return self.properties.metadata(_property)

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
        return self.property.name(pid)

    def property_id(self, name):
        "Obsolete routine kept for compatibility"
        return self.id(name)

    def property_name(self, pid):
        "Obsolete routine kept for compatibility"
        return self.name(pid)

    def property_type(self, _property):
        "Obsolete routine kept for compatibility"
        return self.type(_property)

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
        return self.properties.type(_property)

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
        return self.properties.units(_property)

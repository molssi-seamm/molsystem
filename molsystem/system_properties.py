# -*- coding: utf-8 -*-

"""Property methods for systems."""

import logging
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


class _SystemProperties(object):
    """A class for handling the properties of a system in the star schema."""

    def __init__(self, system):
        self._system = system
        self._sid = system.id
        self._properties = None

    @property
    def system(self):
        """The system these properties are connected with."""
        return self._system

    @property
    def properties(self):
        """The general properties object."""
        if self._properties is None:
            self._properties = self.system.system_db.properties
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

    def get(
        self,
        pattern="*",
        match="glob",
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
        return self.properties.get(
            self._sid,
            is_system=True,
            pattern=pattern,
            match=match,
            include_configuration_properties=include_configuration_properties,
            types=types,
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

    def list(
        self,
        pattern="*",
        match="glob",
        include_configuration_properties=False,
        as_ids=False,
        types=["float", "int", "str", "json"],
    ):
        """Get the list of matching properties for the configuration.

        Parameters
        ----------
        pattern : str = "*"
            The pattern of the property.
        match : str="glob"
            Whether to use exact, glob, or 'like' matching.
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
        return self.properties.list(
            self._sid,
            is_system=True,
            pattern=pattern,
            match=match,
            include_configuration_properties=include_configuration_properties,
            as_ids=as_ids,
            types=types,
        )

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

    def put(self, _property, value):
        """Store the given property value for this system.

        Parameters
        ----------
        _property : int or str
            The id or name of the property.
        value : int, float, or str
            The value to store.
        """
        self.properties.put(self._sid, _property, value, is_system=True)

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

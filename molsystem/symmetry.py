# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)


class _Symmetry(object):
    """A class to handle point and space group symmetry.

    :meta public:
    """

    def __init__(self, configuration):
        """Initialize from the database.

        Parameters
        ----------
        system_db : SystemDB
            The SystemDB instance that we are working with.
        _id : int
            The id of this particular symmetry.
        """
        self._configuration = configuration
        self._system = self._configuration.system
        self._system_db = self._system.system_db
        self._id = configuration.symmetry_id

        super().__init__()

    @property
    def configuration(self):
        """Return the configuration."""
        return self._configuration

    @property
    def cursor(self):
        return self.system_db.cursor

    @property
    def db(self):
        return self.system_db.db

    @property
    def group(self):
        """The point or space group of the system"""
        self.cursor.execute('SELECT "group" FROM symmetry WHERE id = ?', (self.id,))
        return self.cursor.fetchone()[0]

    @group.setter
    def group(self, value):
        self.db.execute(
            'UPDATE symmetry SET "group" = ? WHERE id = ?', (value, self.id)
        )
        self.db.commit()

    @property
    def id(self):
        """The id of this cell."""
        return self._id

    @property
    def system(self):
        """Return the System object that contains this cell."""
        return self._system

    @property
    def system_db(self):
        """Return the SystemDB object that contains this cell."""
        return self._system_db

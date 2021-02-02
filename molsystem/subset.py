# -*- coding: utf-8 -*-
import logging

from .atoms import _SubsetAtoms
from .bonds import _SubsetBonds
from .openbabel import OpenBabelMixin
from .smiles import SMILESMixin
from .template import _Template

logger = logging.getLogger(__name__)


class _Subset(SMILESMixin, OpenBabelMixin, object):
    """:meta public:
    A class providing the API for a subset.

    Parameters
    ----------
    sid : int
        The id of the subset in the subset table.
    logger : logging.Logger = logger
        A logger to use in place of the one from this module.
    """

    def __init__(self, system_db, sid, logger=logger):
        self._system_db = system_db
        self._id = sid
        self._logger = logger

        self._template = None
        self._configuration = None

    @property
    def atoms(self):
        """The atoms for this subset.

        Returns
        -------
        _SubsetAtoms
            The atoms for the subset.
        """
        return _SubsetAtoms(self.configuration, self.id)

    @property
    def bonds(self):
        """The bonds for this subset.

        Returns
        -------
        _SubsetBonds
            The bonds for the subset.
        """
        return _SubsetBonds(self.configuration, self.id)

    @property
    def category(self):
        """The category of this subset."""
        return self.template.category

    @property
    def cursor(self):
        """The a cursor for the database."""
        return self.system_db.cursor

    @property
    def configuration(self):
        """The configuration for this subset."""
        if self._configuration is None:
            sql = "SELECT configuration FROM subset WHERE id = ?"
            self.cursor.execute(sql, (self._id,))
            cid = self.cursor.fetchone()[0]
            self._configuration = self.system_db.get_configuration(cid)
        return self._configuration

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def id(self):
        """Return the id of the subset."""
        return self._id

    @property
    def system_db(self):
        """The system_db that we belong to."""
        return self._system_db

    @property
    def template(self):
        """The template for this subset."""
        if self._template is None:
            sql = "SELECT template FROM subset WHERE id = ?"
            self.cursor.execute(sql, (self._id,))
            tid = self.cursor.fetchone()[0]
            self._template = _Template(self.system_db, tid)
        return self._template

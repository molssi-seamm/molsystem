# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)


class _Symmetry(object):
    """A class to handle point and space group symmetry.

    :meta public:
    """

    def __init__(self, system_db, _id):
        """Initialize from the database.

        Parameters
        ----------
        system_db : SystemDB
            The SystemDB instance that we are working with.
        _id : int
            The id of this particular symmetry.
        """
        self._system_db = system_db
        self._id = _id

        super().__init__()

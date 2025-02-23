# -*- coding: utf-8 -*-

"""A class providing a convenient interface for subsets"""

from itertools import zip_longest
import logging

from .subset import _Subset
from .table import _Table

logger = logging.getLogger(__name__)


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Subsets(_Table):
    """The Subsets class handles the subsets for a single configuration.

    Parameters
    ----------
    configuration : _Configuration
        The configuration to work with.
    logger : Logger
        The logger to use, defaults to the one for this module.
    """

    def __init__(self, configuration, logger=logger):
        self._configuration = configuration
        self.logger = logger

        self._cid = configuration.id
        self._system_db = configuration.system_db

        super().__init__(self._system_db, "subset")

        self._sa_table = self._system_db["subset_atom"]

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    @property
    def n_subsets(self):
        """The number of subsets for a configuration.

        Returns
        -------
        int
            The number of subsets in the configuration.
        """
        sql = """\
        SELECT COUNT(*)
          FROM subset
         WHERE configuration = ?
        """

        self.cursor.execute(sql, (self._cid,))
        result = self.cursor.fetchone()[0]
        return result

    def append(self, n=None, **kwargs):
        """Append one or more rows

        The keywords are the names of attributes and the value to use.
        The default value for any attributes not given is used unless it is
        'None' in which case an error is thrown. It is an error if there is not
        an exisiting attribute corresponding to any given as arguments.

        Parameters
        ----------
        n : int = None
            The number of rows to create, defaults to the number of items
            in the longest attribute given.
        kwargs :
            any number <attribute name> = <value> keyword arguments giving
            existing attributes and values.

        Returns
        -------
        [int]
            The ids of the created rows.
        """
        kwargs["configuration"] = self._cid
        ids = super().append(n, **kwargs)

        return ids

    def create(self, template, atoms=None, templateatoms=None):
        """Create a subset given a template and optionally atoms.

        Parameters
        ----------
        template : int or _Template
            The template for the subset
        atoms : [int] = None
            Optional list of atom ids to connect to the subset.
        templateatoms : [int] = None
            Optional list of template atoms to connect to the atoms
            in the subset.

        Returns
        -------
        _Subset
            The subset.
        """
        if isinstance(template, int):
            tid = template
        else:
            tid = template.id
        sid = self.append(template=tid)[0]

        if atoms is not None:
            if templateatoms is not None:
                self._sa_table.append(
                    subset=sid, atom=atoms, templateatom=templateatoms
                )
            else:
                self._sa_table.append(subset=sid, atom=atoms)

        return _Subset(self.system_db, sid)

    def delete(self, ids):
        """Remove one or more subsets

        Parameters
        ----------
        ids : [int] of 'all'
            The subsets to delete.

        Returns
        -------
        None
        """
        if ids == "all":
            sql = """\
            DELETE FROM subset
             WHERE configuration = ?
            """
            self.db.execute(sql, (self._cid,))
        else:
            if isinstance(ids, int):
                self.db.execute("DELETE FROM subset WHERE id = ?", (ids,))
            else:
                self.db.executemany("DELETE FROM subset WHERE id = ?", ids)

    def generate(self, template):
        """Generate the subsets matching a full template.

        Parameters
        ----------
        template : int or _Template
           The template.

        Returns
        -------
        [_Subset]
            The subsets.
        """
        pass

    def get(self, template):
        """Get the subsets given a template.

        Parameters
        ----------
        template : int or _Template
           The template.

        Returns
        -------
        [_Subset]
            The subsets.
        """
        if isinstance(template, int):
            tid = template
        else:
            tid = template.id

        return [_Subset(self.system_db, x) for x in self.get_ids(tid)]

    def get_ids(self, template):
        """Find ids of subsets given a template.

        Parameters
        ----------
        template : int or _Template
            The template for the subset

        Returns
        -------
        [int]
            The ids of the subsets, and empty list if there are none.
        """
        if isinstance(template, int):
            tid = template
        else:
            tid = template.id

        sql = """
        SELECT id FROM subset
         WHERE template = ?
           AND configuration = ?
        """
        return [x[0] for x in self.db.execute(sql, (tid, self._cid))]

    def get_counts(self, configuration="all"):
        """Get the counts of subsets per template for configurations.

        Parameters
        ----------
        configuration : int or MolSystem._Configuration or "all"
            The configuration, as an object or its id, or "all" for all
            configurations. Default is "all".

        Returns
        -------
        [[configuration_id, name, count]]
        """

        if configuration == "all":
            sql = (
                "SELECT subset.configuration, name, count(1) FROM subset, template "
                " WHERE template = template.id "
                " GROUP BY subset.configuration, name "
                " ORDER BY subset.configuration"
            )
            return [[*x] for x in self.db.execute(sql)]
        else:
            raise NotImplementedError("get_counts only take 'all' configurations")

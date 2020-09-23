# -*- coding: utf-8 -*-

"""A class providing a convenient interface for subsets
"""

from itertools import zip_longest
import logging
from typing import TypeVar, Dict, Any

from molsystem.table import _Table as Table

System_tp = TypeVar("System_tp", "System", None)
Templates_tp = TypeVar("Templates_tp", "_Templates", str, None)

logger = logging.getLogger(__name__)


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,...s3n-1), ..."
    return zip_longest(*[iter(iterable)] * n)


class _Subsets(Table):
    """The Subset class works with the tables controlling subsets.

    See the main documentation of SEAMM for a detailed description of the
    database scheme underlying the system and hence the subsets.. The
    following tables handle subsets:

    """

    def __init__(self, system: System_tp, tablename: str = 'subset') -> None:
        super().__init__(system, tablename)

        self._configuration_subset_table = self.system['configuration_subset']

    def n_subsets(self, configuration=None):
        """The number of subsets for a configuration.

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        int
            The number of subsets in the configuration.
        """
        if configuration is None:
            configuration = self.system.current_configuration

        self.cursor.execute(
            f'SELECT COUNT(*) FROM {self.table}, "configuration_subset"'
            '  WHERE id = subset AND configuration = ?', (configuration,)
        )
        result = self.cursor.fetchone()[0]
        return result

    def append(
        self,
        n: int = None,
        configuration: int = None,
        **kwargs: Dict[str, Any]
    ) -> None:
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
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.
        kwargs :
            any number <attribute name> = <value> keyword arguments giving
            existing attributes and values.

        Returns
        -------
        [int]
            The ids of the created rows.
        """
        if configuration is None:
            configuration = self.system.current_configuration

        ids = super().append(n, **kwargs)

        # and link to the configuration
        self._configuration_subset_table.append(
            configuration=configuration, subset=ids
        )

        return ids

    def delete(self, ids, configuration=None):
        """Remove one or more subsets

        Parameters
        ----------
        ids : [int]
            The subsets to delete.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        None
        """
        if ids == 'all':
            if configuration is None:
                configuration = self.system.current_configuration
            sql = (
                f'DELETE FROM {self.table}, "configuration_subset"'
                '  WHERE id = subset and configuration = ?'
            )
            self.db.execute(sql, (configuration,))
        else:
            if isinstance(ids, int):
                self.db.execute(
                    f"DELETE FROM {self.table} WHERE id = ?", (ids,)
                )
            else:
                self.db.executemany(
                    f"DELETE FROM {self.table} WHERE id = ?", ids
                )

    def create(
        self, template, configuration=None, atoms=None, templateatoms=None
    ):
        """Create a subset given a template and optionally atoms.

        Parameters
        ----------
        template : int
            The template for the subset
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.
        atoms : [int] = None
            Optional list of atom ids to connect to the subset.
        templateatoms : [int] = None
            Optional list of template atoms to connect to the atoms
            in the subset.

        Returns
        -------
        int
            The id of the subset.
        """
        sid = self.append(template=template, configuration=configuration)[0]

        if atoms is not None:
            sa = self.system['subset_atom']
            if templateatoms is not None:
                sa.append(subset=sid, atom=atoms, templateatom=templateatoms)
            else:
                sa.append(subset=sid, atom=atoms)

        return sid

    def find(self, template, configuration=None):
        """Find subsets given a template.

        Parameters
        ----------
        template : int
            The template for the subset
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        [int]
            The ids of the subsets, and empty list if there are none.
        """
        if configuration is None:
            configuration = self.system.current_configuration

        result = []
        sql = (
            f'SELECT id FROM {self.table}, "configuration_subset"'
            '  WHERE id = subset AND configuration = ? AND template = ?'
        )
        for row in self.db.execute(sql, (configuration, template)):
            result.append(row['id'])

        return result

    def template(self, sid):
        """The template for the given subset.

        Parameters
        ----------
        sid : int
            The id of the subset.

        Returns
        -------
        int
            The id of the associated template
        """
        self.cursor.execute(
            f'SELECT template FROM {self.table} WHERE id = ?', (sid,)
        )
        return self.cursor.fetchone()[0]

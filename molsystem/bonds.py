# -*- coding: utf-8 -*-

"""A dictionary-like object for holding bonds

Based on tables in an SQLite database.
"""

import logging
import sqlite3
from typing import Any, Dict, TypeVar

import numpy as np
import pandas

from molsystem.table import _Table as Table
from molsystem.frozencolumn import _FrozenColumn as FrozenColumn
from molsystem.column import _Column as Column

System_tp = TypeVar("System_tp", "System", "Bonds", None)

logger = logging.getLogger(__name__)


class _Bonds(Table):
    """The Bonds class describes the bonds in the system.

    Three attributes are required: the atoms 'i' and 'j' of the bond, and the
    bond order, which defaults to a single bond.

    Since there are two possible representations for a bond, i-j and j-i,
    the bonds are stored with i < j in order to make the indexing unique.

    Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Bonds can be added ('append') or removed ('delete').

    In order to handle changes in bonding (and numbers of atoms) the
    system class has several tables covering the atoms and bonds. See
    the main documentation for a more detailed description. The key
    tables are:

    configuration -- a list of configurations in e.g. a trajectory
    configuration_subset -- connnects configurations with subset

    template -- a simple list of templates
    templateatom -- holds the atoms in a template, if present.
    templatebond -- holds the bonds between the template atoms, if present

    subset -- an instantiation of a template, connected with one or more
              configurations of the system.
    subset_atom -- connects the subset to the atoms, and optionally connects
                   the atom to a template atom.

    atom -- the nonvarying part of the description of the atoms
    coordinates -- the varying part of the description of atoms

    This class handles the bonds for the main system, but it does so through
    subset "all", which is a special subset that contains all of the atoms in
    the system. This subset is defined by one or more templates "all" which are
    instantiated in the subset table and connected to configurations via the
    configuration_subset table.

    Underneath, this is handled by a subset "all", which is a special
    subset that contains all of the atoms in the system. This subset is
    defined by one or more templates "all" which are instantiated in the
    subset table and connected to configurations via the
    configuration_subset table.

    If the number of atoms in the system is fixed and there are no bonds
    or they don't change (the case in classical MD), there will be one
    template "all" instantiated as one subset. All configurations will
    point to this one subset. If there are no bonds in the system, then
    this is all that is needed. However if there are bonds, the atoms
    involved in bonds need to be in the templateatom table and the bonds
    entered in the templatebond table.

    If the bonds change over time, then for each bond configuration
    there needs to be a distinct "all" template, with the appropriate
    templateatoms and templatebonds. This template is instantiated as a
    new subset and connected to one or more configurations.

    If the number of atoms changes over time, each new set of atoms also
    requires a new template and subset connected to the atoms in that
    configuration. A subset may persist over a number of configurations,
    but each time the atoms change a new subset is required.

    If both the atoms and the bonds change, then a corresponding subset
    is needed each time either changes.
    """

    def __init__(self, system: System_tp, table: str = 'templatebond') -> None:

        super().__init__(system, table)

        self._bond_db = None

    def __getitem__(self, key):
        return self.get_column(key)

    def __eq__(self, other) -> Any:
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    @property
    def bond_db(self):
        if self._bond_db is None:
            self._bond_db = sqlite3.connect(self._system._filename)
            self._bond_db.row_factory = self._row_factory
            self._bond_db.execute('PRAGMA foreign_keys = ON')
        return self._bond_db

    def append(
        self,
        configuration: int = None,
        bonds=None,
        **kwargs: Dict[str, Any]
    ) -> None:
        """Append one or more bonds

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.

        This is where it becomes a bit complicated. The atoms i & j are atoms
        in the main 'atom' table; however, the bonds are stored in
        'templatebond' referencing atoms in 'templateatom' for the template
        'all'. The real atoms are connected to the template atoms via the
        'subset_atom' table, which has an extra field referencing the template
        atom.
        """

        if bonds is not None:
            if len(kwargs) > 0:
                raise RuntimeError(
                    'Please do not give both bonds and arrays of data!'
                )
            if isinstance(bonds, sqlite3.Row):
                # One bond
                kwargs = {}
                for key, value in zip(bonds.keys(), bonds):
                    kwargs[key] = value
            else:
                first = True
                for bond in bonds:
                    if first:
                        keys = bond.keys()
                        for key, value in zip(keys, bond):
                            kwargs[key] = [value]
                        first = False
                    else:
                        for key, value in zip(keys, bond):
                            kwargs[key].append(value)

        # Check keys and lengths of added bonds
        if 'i' not in kwargs or 'j' not in kwargs:
            raise KeyError("The atoms i & j are required!")

        i = kwargs.pop('i')
        j = kwargs.pop('j')
        try:
            len_i = len(i)
        except TypeError:
            len_i = 1
            i = (i,)
        try:
            len_j = len(j)
        except TypeError:
            len_j = 1
            j = (j,)

        if len_i == 0 and len_j == 0:
            return

        if len_i == 1:
            if len_j > 1:
                i = len_j * [i[0]]
        elif len_j == 1:
            j = len_i * [j[0]]
        elif len_i != len_j:
            raise IndexError(
                f'key "j" has the wrong number of values, {len_j}. '
                f'Should be 1 or the number of values in i, {len_i}.'
            )
        # Ensure that i < j

        if isinstance(i, np.ndarray):
            i = i.tolist()
        if isinstance(j, np.ndarray):
            j = j.tolist()

        i2 = []
        j2 = []
        for i_, j_ in zip(i, j):
            # will need to handle offsets here at some point
            if not isinstance(i_, int) or not isinstance(j_, int):
                raise TypeError(
                    "'i' and 'j', the atom indices, must be integers"
                )
            if i_ < j_:
                i2.append(i_)
                j2.append(j_)
            else:
                i2.append(j_)
                j2.append(i_)

        # Get the template atoms corresponding to the i and j atoms,
        # adding them if needed.
        if configuration is None:
            configuration = self._system.current_configuration
        subset = self._system.all_subset(configuration)
        template = self._system.all_template(configuration)

        # First get current map of atoms to templateatoms
        map = {}
        for row in self.db.execute(
            "SELECT atom, templateatom FROM subset_atom"
            " WHERE subset = ?", (subset,)
        ):
            templateatom = row['templateatom']
            if templateatom is not None:
                map[row['atom']] = templateatom

        # Use dictionary because there may be duplicates
        missing = {}
        for i_ in i2:
            if i_ not in map:
                missing[i_] = None
        for j_ in j2:
            if j_ not in map:
                missing[j_] = None

        if len(missing) > 0:
            # Need to add to the all template.
            table = self._system['templateatom']
            templateatoms = table.append(n=len(missing), template=template)

            parameters = []
            # The sorting here is to keep the atoms in the same order,
            # which is what the user expects.
            for i_, tatom in zip(sorted(missing), templateatoms):
                map[i_] = tatom
                parameters.append((tatom, i_, subset))
            self.cursor.executemany(
                "UPDATE subset_atom SET templateatom = ?"
                " WHERE atom = ? AND subset = ?", parameters
            )

        # get the lists of template atoms for the bonds
        ti = []
        tj = []
        for i_ in i2:
            ti.append(map[i_])
        for j_ in j2:
            tj.append(map[j_])

        # and ... finally ... add the bonds
        if 'bondorder' in kwargs:
            bondorders = kwargs.pop('bondorder')
            self._system['templatebond'].append(
                i=ti, j=tj, bondorder=bondorders
            )
        else:
            self._system['templatebond'].append(i=ti, j=tj)

    def bonds(self, subset=None, configuration=None):
        """Returns an iterator over the rows of the bonds.

        Note that these are references to the main atom table.

        Parameters
        ----------
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        sqlite3.Cursor
            A cursor that returns sqlite3.Row objects for the bonds.
        """
        if subset is None:
            subset = self.system.all_subset(configuration)

        columns = []
        for column in self.attributes:
            if column != 'i' and column != 'j':
                columns.append(f'templatebond."{column}"')
        column_defs = ', '.join(columns)
        sql = (
            f"SELECT iatom.atom as i, jatom.atom as j, {column_defs}"
            "   FROM templatebond, subset_atom as iatom, subset_atom as jatom"
            "  WHERE templatebond.i = iatom.templateatom"
            "    AND templatebond.j = jatom.templateatom"
            "    AND iatom.subset = ? and jatom.subset = ?"
        )
        return self.db.execute(sql, (subset, subset))

    def contains_bond(self, key):
        if isinstance(key, self.bond_tuple):
            i = key.i
            j = key.j
        else:
            i, j = key
        # get canonical order
        if i > j:
            j, i = i, j
        sql = (
            "SELECT COUNT(*) FROM templatebond, subset_atom as iatom,"
            "              subset_atom as jatom"
            " WHERE iatom.subset = ? AND jatom.subset = ?"
            "   AND templatebond.i = iatom.templateatom"
            "   AND templatebond.j = jatom.templateatom"
            "   AND iatom.atom = ? AND jatom.atom = ?"
        )
        all_subset = self._system.all_subset()
        self.cursor.execute(sql, (all_subset, all_subset, i, j))
        return self.cursor.fetchone()[0] == 1

    def delete_bond(self, i, j):
        """Remove a bond, if present"""
        # get canonical order
        if i > j:
            j, i = i, j
        sql = (
            "SELECT templatebond.i as i, templatebond.j as j"
            "  FROM templatebond, subset_atom as iatom, subset_atom as jatom"
            " WHERE iatom.subset = ? AND jatom.subset = ?"
            "   AND templatebond.i = iatom.templateatom"
            "   AND templatebond.j = jatom.templateatom"
            "   AND iatom.atom = ? AND jatom.atom = ?"
        )
        all_subset = self._system.all_subset()
        self.cursor.execute(sql, (all_subset, all_subset, i, j))
        row = self.cursor.fetchone()
        if row is not None:
            sql = "DELETE FROM templatebond WHERE i = ? AND j = ?"
            self.cursor.execute(sql, row)

    def get_bond(self, i, j):
        # get canonical order
        if i > j:
            j, i = i, j
        sql = (
            "SELECT templatebond.* FROM templatebond, subset_atom as iatom,"
            "              subset_atom as jatom"
            " WHERE iatom.subset = ? AND jatom.subset = ?"
            "   AND templatebond.i = iatom.templateatom"
            "   AND templatebond.j = jatom.templateatom"
            "   AND iatom.atom = ? AND jatom.atom = ?"
        )
        all_subset = self._system.all_subset()
        cursor = self.db.execute(sql, (all_subset, all_subset, i, j))
        row = cursor.fetchone()
        if row is None:
            raise KeyError(f'No bond from {i} to {j} found')
        return row

    def get_column(self, key, configuration=None):
        """Return a column from the templatebonds table.

        Parameters
        ----------
        key : str
            The column (attribute) to getLogger
        configuration : int = None
            Get the bonds for this configuration, defaults to the current
            configuration.

        Returns
        -------
        column : Column
            The column object requested.
        """
        if key == 'i' or key == 'j':
            all_subset = self._system.all_subset(configuration)
            if key == 'i':
                col = 'iatom.atom'
            else:
                col = 'jatom.atom'
            sql = (
                f"SELECT {col}"
                "  FROM templatebond,"
                "       subset_atom as iatom,"
                "       subset_atom as jatom"
                " WHERE templatebond.i = iatom.templateatom"
                "   AND templatebond.j = jatom.templateatom"
                f"  AND iatom.subset = {all_subset}"
                f"  AND jatom.subset = {all_subset}"
            )
            table = Table(self._system, 'atom')
            return FrozenColumn(table, 'id', sql)
        else:
            all_template = self._system.all_template(configuration)
            sql = (
                f'SELECT templatebond.rowid, templatebond."{key}"'
                "  FROM templatebond,"
                "       templateatom as iatom,"
                "       templateatom as jatom"
                " WHERE templatebond.i = iatom.id"
                "   AND templatebond.j = jatom.id"
                f"  AND iatom.template = {all_template}"
                f"  AND jatom.template = {all_template}"
            )
            table = Table(self._system, 'templatebond')
            return Column(table, key, sql=sql)

    def n_bonds(self, subset: int = None, configuration: int = None) -> int:
        """The number of bonds.

        Parameters
        ----------
        bonds : int = None
            Get the bonds for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        int
            Number of bonds
        """
        if subset is None:
            template = self.system.all_template(configuration)
        else:
            template = self.system.subsets.template(subset)

        self.cursor.execute(
            "SELECT COUNT(*) FROM templateatom, templatebond "
            " WHERE templatebond.i = templateatom.id"
            "   AND templateatom.template = ?", (template,)
        )
        return self.cursor.fetchone()[0]

    def remove(self, atoms=None, subset=None, configuration=None):
        """Removes all the bonds for the atoms, or in a subset or configuration

        Parameters
        ----------
        atoms : [int] = None
            The list of atoms whose bonds are removed
        subset : int = None
            Get the atoms for the subset. Defaults to the 'all/all' subset
            for the configuration given.
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        None
        """
        if atoms is not None:
            if subset is None:
                subset = self.system.all_subset(configuration)
            # Delete the bonds in the template
            sql = (
                "DELETE FROM templatebond"
                " WHERE i in ("
                "     SELECT templateatom FROM subset_atom"
                "      WHERE subset = ? AND atom = ?"
                " ) OR j in ("
                "     SELECT templateatom FROM subset_atom"
                "      WHERE subset = ? AND atom = ?"
                " )"
            )
            parameters = [(subset, i, subset, i) for i in atoms]
            self.db.executemany(sql, parameters)

            # and the atoms...
            sql = (
                "DELETE FROM templateatom"
                " WHERE id in ("
                "     SELECT templateatom FROM subset_atom"
                "      WHERE subset = ? AND atom = ?"
                " )"
            )
            parameters = [(subset, i) for i in atoms]
            self.db.executemany(sql, parameters)
        else:
            if subset is None:
                subset = self.system.all_subset(configuration)

            sql = (
                "DELETE FROM templatebond"
                " WHERE i in ("
                "     SELECT id FROM templateatom, subset_atom"
                "      WHERE id = templateatom AND subset = ?"
                " ) AND j in ("
                "     SELECT id FROM templateatom, subset_atom"
                "      WHERE id = templateatom AND subset = ?"
                " )"
            )
            self.db.execute(sql, (subset, subset))

    def to_dataframe(self, configuration=None):
        """Return the bonds as a Pandas Dataframe."""
        all_subset = self._system.all_subset(configuration)
        columns = []
        for column in self.attributes:
            if column != 'i' and column != 'j':
                columns.append(f'templatebond."{column}"')
        column_defs = ', '.join(columns)
        sql = (
            "SELECT templatebond.rowid,"
            "       iatom.atom as i,"
            "       jatom.atom as j,"
            f"      {column_defs}"
            "  FROM templatebond, subset_atom as iatom, subset_atom as jatom"
            " WHERE templatebond.i = iatom.templateatom"
            "   AND templatebond.j = jatom.templateatom"
            "   AND iatom.subset = ? and jatom.subset = ?"
        )

        data = {}
        for line in self.cursor.execute(sql, (all_subset, all_subset)):
            data[line[0]] = line[1:]
        columns = [x[0] for x in self.cursor.description[1:]]

        df = pandas.DataFrame.from_dict(data, orient='index', columns=columns)

        return df

        # columns = [*self.attributes]
        # df = pandas.DataFrame.from_records(
        #     self.bonds(configuration), columns=columns
        # )
        # return df

    def _row_factory(self, cursor, row):
        return self.bond_tuple(*row)

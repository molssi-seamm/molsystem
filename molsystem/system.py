# -*- coding: utf-8 -*-

"""A dictionary-like object for holding a system
"""

from collections.abc import MutableMapping
import logging
import pprint  # noqa: F401

from .cif import CIFMixin
from .configuration import _Configuration
from .table import _Table

logger = logging.getLogger(__name__)


class _System(CIFMixin, MutableMapping):
    """A single system -- molecule, crystal, etc. -- in SEAMM.

    Based on a SQLite database, this class provides a general
    datastructure for storing a 'system', which is an atom, molecule,
    crystal or a collection of any number of such components all held in
    a single system, in analogy to a thermodynamic system, which is what
    an experiment is perfromed on.

    The system can hold multiple different configurations of the atoms,
    such as found in an MD trajectory or conformer search. The number of
    atoms may change across configurations as may the bonding, so the
    system is capabale of handling a trajectory for e.g. a grand
    canonical simulation or simulating desposition on a surface. The
    system also can contain bonds as are used in classical valence
    forcefields, and can also cope with changing bonds in e.g. reactive
    forcefields.

    The system class also supports the concept of template molecules or
    subsets, which can be used to define molecules, residues, etc. which
    are then repeated in the full structure. Templates make the process
    of building complex, repeated structures simpler and more efficient,
    and also help track the individual components. There is a lighter
    variation of templates called subsets which do not define a
    substructure, but can label and track such substructures. Atoms may
    belong to more than one subset, or to none other than the special
    subset 'all' which contains all the atoms in a given configuration.

    The key tables are:

    system -- a container for all the components, plus storage for the
        system-wide parameters and options.

    configuration -- a list of configurations in e.g. a trajectory.
    configuration_subset -- connnects configurations with subsets.

    cell -- the information about the periodicity of a configuration, if
        needed.

    symmetry -- the point or space group symmetry of the configuration.
    symop -- the symmetry operations for the group.

    template -- a list of all templates.

    subset -- an instantiation of a template, connected with one or more
        configurations of the system.
    subset_atom -- connects the subset to the atoms, and optionally connects
        the atom to a template atom.

    atom -- the nonvarying part of the description of the atoms
    coordinates -- the varying part of the description of atoms

    element -- basic information about the elements. The atomic number
        (atno) is used for foreign keys in other tables.

    This system class provides an anchor for the other component and
    handles system-wide information such as the periodicity and the
    coordinate system used.

    It also tracks the current configuration and template, which are
    conveniences to avoid having to specify the configuration or
    template for every access to the data.

    There is one special, required subset 'all' that contains all of the
    atoms in each configuration of the system. If there are differing
    numbers or types of atoms in the system, or if the bonding changes,
    there will be multiple versions of the subset, each corresponding to
    a different set of atoms or bonds.

    One or more configurations are connected with each 'all' subset,
    which is how the atoms and bonding are connected to the
    configurations.

    :meta public:
    """

    def __init__(self, system_db, _id, logger=logger):
        self._system_db = system_db
        self._id = _id
        self.logger = logger
        self._current_configuration_id = None  # The current configuration
        self._checkpoints = []
        self._items = {}

    def __del__(self):
        """Delete all instance variables that are objects."""
        del self._system_db

    def __getitem__(self, key):
        """Allow [] access to the dictionary!

        Because some of the items, such as template contain state,
        we need to ensure that the same object is used everywhere. Hence
        the self._items array to store the instances.

        Parameters
        ----------
        key : str
            The table name or name of a multi-table item like atoms.

        Returns
        -------
        table : Table or Atom, Template, Templateatoms, ...
        """
        if key not in self._items:
            self._items[key] = _Table(self, key)
        return self._items[key]

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        raise NotImplementedError(f"Table '{key}' cannot be created yet")

    def __delitem__(self, key):
        """Allow deletion of keys"""
        if key in self:
            self.cursor.execute(f"DROP TABLE '{key}'")

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self.list())

    def __len__(self):
        """The len() command"""
        self.cursor.execute(
            "SELECT COUNT(*)" "  FROM sqlite_master" " WHERE type = 'table'"
        )
        return self.cursor.fetchone()[0]

    # def __repr__(self):
    #     """The string representation of this object"""
    #     raise NotImplementedError()

    # def __str__(self):
    #     """The pretty string representation of this object"""
    #     raise NotImplementedError()

    def __contains__(self, item):
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.

        tmp_item = item.strip('"')
        self.cursor.execute(
            "SELECT COUNT(*)"
            "  FROM sqlite_master"
            f" WHERE type = 'table' and name = '{tmp_item}'"
        )
        return self.cursor.fetchone()[0] == 1

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        if self is other:
            return True

        tables = set(self.list())
        other_tables = set(other.list())

        added = tables - other_tables
        deleted = other_tables - tables
        in_common = tables & other_tables

        if len(added) > 0:
            return False
        if len(deleted) > 0:
            return False

        # Need the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        name = self.name
        other_name = other.name
        detach = False
        if name != other_name and not self.is_attached(other_name):
            # Attach the other system in order to do comparisons.
            self.attach(other)
            detach = True

        # Check the tables in both systems
        result = True
        for table in in_common:
            if _Table(self, table) != _Table(other, table):
                result = False
                break

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    @property
    def configuration(self):
        """The configuration object for the current configuration."""
        if self.n_configurations == 0:
            return None
        if self._current_configuration_id is None:
            self._current_configuration_id = self.configuration_ids[-1]
        return _Configuration(
            _id=self._current_configuration_id, system_db=self.system_db
        )

    @configuration.setter
    def configuration(self, value):
        if isinstance(value, _Configuration):
            value = value.id
        if value not in self.configuration_ids:
            raise KeyError(f"configuration '{value}' does not exist.")
        self._current_configuration_id = value

    @property
    def configuration_names(self):
        """Get the names of the configurations.

        Returns
        -------
        [str]
            The names of the configurations.
        """
        # The name of the configuration ... find it.
        sql = "SELECT name FROM configuration WHERE system = ?"
        return [x[0] for x in self.db.execute(sql, (self.id,))]

    @property
    def configurations(self):
        """The list of configuration objects for this system."""
        result = []
        for _id in self.configuration_ids:
            result.append(_Configuration(system_db=self.system_db, _id=_id))
        return result

    @property
    def configuration_ids(self):
        """The list of configuration ids."""
        result = []
        sql = "SELECT id FROM configuration WHERE system = ?"
        for row in self.db.execute(sql, (self.id,)):
            result.append(row[0])
        return result

    @property
    def cursor(self):
        return self.system_db.cursor

    @property
    def db(self):
        return self.system_db.db

    @property
    def id(self):
        """The id of this system."""
        return self._id

    @property
    def name(self):
        """Return the name of this system."""
        self.cursor.execute("SELECT name FROM system WHERE id = ?", (self.id,))
        result = self.cursor.fetchone()
        if result is None:
            return None
        else:
            return result[0]

    @name.setter
    def name(self, value):
        self.db.execute("UPDATE system SET name = ? WHERE id = ?", (value, self.id))
        self.db.commit()

    @property
    def n_configurations(self):
        """The number of configurations of the system."""
        sql = "SELECT COUNT(*) FROM configuration WHERE system = ?"
        self.cursor.execute(sql, (self.id,))
        return self.cursor.fetchone()[0]

    @property
    def system_db(self):
        """Return the SystemDB object that contains this system."""
        return self._system_db

    def copy_configuration(
        self,
        configuration=None,
        name=None,
    ):
        """Add a new configuration by copying another configuration.

        Parameters
        ----------
        configuration : int = None
            The configuration to copy. Defaults to the current configuration.
        name : str = None
            A textual name for the configuration (optional)

        Returns
        -------
        cid : int
            The id of the new configuration.

        """
        configuration = self.get_configuration(configuration)

        cid = self["configuration"].append(
            system=self.id,
            name=name,
            periodicity=configuration.periodicity,
            coordinatesystem=configuration.coordinatesystem,
            symmetry=configuration.symmetry_id,
            cell=configuration.cell_id,
            atomset=configuration.atomset,
            bondset=configuration.bondset,
        )[0]

        return cid

    def create_configuration(
        self,
        name="",
        periodicity=0,
        coordinatesystem=None,
        symmetry=None,
        cell_id=None,
        atomset=None,
        bondset=None,
        make_current=True,
    ):
        """Add a new configuration to the system.

        Parameters
        ----------
        name : str = None
            A textual name for the configuration (optional)
        periodicity : int = 0
            The periodicity, 0, or 3 for 0-D or 3-D at the moment
        coordinatesystem : str = None
            The coordinate system, 'Cartesian' or 'fractional', to use.
            Defaults to Cartesian for molecules and fractional for crystals.
        symmetry : int or str = None
            The id or name of the point or space group (optional)
        cell_id : int = None
            The id of the _Cell object
        atomset : int = None
            The set of atoms in this configurationn
        bondset : int = None
            The bonds in this configuration
        make_current : bool = True
            If True, make the current configuration.

        Returns
        -------
        _Configuration
            The new configuration.
        """
        kwargs = {}
        if name is not None:
            kwargs["name"] = name
        if periodicity != 0:
            if coordinatesystem is None:
                kwargs["coordinatesystem"] = "fractional"
            else:
                if coordinatesystem.lower()[0] == "c":
                    kwargs["coordinatesystem"] = "Cartesian"
                else:
                    kwargs["coordinatesystem"] = "fractional"
        if symmetry is not None:
            kwargs["symmetry"] = symmetry
        if cell_id is not None:
            kwargs["cell"] = cell_id
        if atomset is not None:
            kwargs["atomset"] = atomset
        if bondset is not None:
            kwargs["bondset"] = bondset

        cid = self["configuration"].append(system=self.id, **kwargs)[0]
        configuration = _Configuration(_id=cid, system_db=self.system_db)

        # If the atomset was given, need to create dummy coordinates
        if atomset is not None:
            n_atoms = configuration.n_atoms
            data = {"configuration": configuration.id}
            data["atom"] = configuration.atoms.ids
            data["x"] = [0.0] * n_atoms
            data["y"] = [0.0] * n_atoms
            data["z"] = [0.0] * n_atoms
            table = _Table(self.system_db, "coordinates")
            table.append(n=n_atoms, **data)

        if make_current:
            self._current_configuration_id = cid

        return configuration

    def create_table(self, name, cls=_Table, other=None):
        """Create a new table with the given name.

        Parameters
        ----------
        name : str
            The name of the new table.

        cls : _Table subclass
            The class of the new table, defaults to _Table

        Returns
        -------
        table : class _Table
            The new table
        """
        if name in self:
            raise KeyError(f"'{name}' already exists in the system.")

        self._items[name] = cls(self, name, other)
        return self._items[name]

    def diff(self, other):
        """Differences between this system and another."""
        result = {}

        if self is other:
            return result

        tables = set(self.list())
        other_tables = set(other.list())

        added = tables - other_tables
        deleted = other_tables - tables
        in_common = tables & other_tables

        if len(added) > 0:
            result["tables added"] = added
        if len(deleted) > 0:
            result["tables deleted"] = deleted

        # Need the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        name = self.nickname
        other_name = other.nickname
        detach = False
        if name != other_name and not self.is_attached(other_name):
            # Attach the other system in order to do comparisons.
            self.attach(other)
            detach = True

        # Check the tables in both systems
        for table in in_common:
            table1 = _Table(self, table)
            table2 = _Table(other, table)
            tmp = table1.diff(table2)
            if len(tmp) > 0:
                result[f"table '{table}' diffs"] = tmp

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    def get_configuration(self, cid):
        """Get the specified configuration object.

        Parameters
        ----------
        cid : int or str
            The id (int) or name (str) of the configuration

        Returns
        -------
        _Configuration
            The requested configuration.

        Raises
        ------
        ValueError
            If the configuration does not exist, or more than one have the
            requested name.
        """
        if isinstance(cid, str):
            # The name of the configuration ... find it.
            cid = self.get_configuration_id(cid)

        return _Configuration(_id=cid, system_db=self.system_db)

    def get_configuration_id(self, name):
        """Get the id of the specified configuration.

        Parameters
        ----------
        name : str
            The name of the configuration

        Returns
        -------
        int
            The id of the requested configuration.

        Raises
        ------
        ValueError
            If the configuration does not exist, or more than one have the
            requested name.
        """
        # The name of the configuration ... find it.
        sql = "SELECT id FROM configuration WHERE system = ? AND name = ?"
        ids = [x[0] for x in self.db.execute(sql, (self.id, name))]
        if len(ids) == 0:
            raise ValueError(f"The configuration '{name}' does not exist.")
        elif len(ids) > 1:
            raise ValueError(f"There is more than one configuration named '{name}'")
        return ids[0]

    def list(self):
        """Return a list of all the tables in the system."""
        result = []
        for row in self.db.execute(
            "SELECT name" "  FROM sqlite_master" " WHERE type = 'table'"
        ):
            result.append(row["name"])
        return result

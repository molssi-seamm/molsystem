# -*- coding: utf-8 -*-

"""A dictionary-like object for holding a system
"""

import collections.abc
import logging
import sqlite3

import molsystem
from .cif import CIFMixin
from .configuration import _Configuration
from .properties import _Properties
from .system import _System
from .table import _Table
from .templates import _Templates

logger = logging.getLogger(__name__)


class SystemDB(CIFMixin, collections.abc.MutableMapping):
    """A database of systems for SEAMM.

    A class based on a SQLite database for describing molecular and
    periodic systems.

    See the documentation at https://molssi-seamm.github.io/molsystem
    for a complete description of the database schema, including a
    diagram.

    The key concepts are:

    System
        The overall container for the configurations, atoms, bonds,
        etc. A system has one or more configurations (conformers),
        each of which details the atoms, bonds, etc.

    Configuration
        The configuration is a single instance of the system, with
        atoms, coordinates, bonds and subsets. It is close to what
        most programs consider the molecule or crystal.

        There may be different set of atoms in different
        configurations of a single system. This supports e.g. grand
        canonical ensembles.

        The set of bonds may also differ between different
        configurations, whether or not the set of atoms changes. This
        supports e.g. reactive forcefields.

        Each configuration also has its own set of subsets, which are
        a way of collecting atoms into groups, and can be thought of
        as a generalization of the chain and residue nomenclature for
        proteins. A configuration may have any number of subsets, and
        atoms may be in more than one subset.

    Atoms
        Each configuration contains a set of atoms, though different
        configurations may contain different atoms. Each atom is
        identified by its atomic number and unique name and has
        coordinates. It may also have other attributes, but that
        depends on the simulation.

    Bonds
        Each configuration may also have information on bonds between
        the atoms. A bond connects two atoms, and has a bond order
        (single, double, triple, aromatic, ...). In periodic systems
        with an infinite network of bonds, each bond also has a cell
        offset to identify the relative cells of the two atoms.

    Templates and Subsets
        Subsets define groups of atoms in a general way. They are
        defined by a template, which may be nothing more than a type
        and name, e.g. 'residue/ala', or the template may be linked
        with a system which has atoms, bonds, etc.

        Each subset is linked to its group atoms. For example, the
        template 'residue/ala' mentioned above would have a subset for
        each alanine in the protein. If the template were linked to a
        system which was the alanine residue, then each atom would
        also be connected with the appropriate atom in the template
        system, so even if the order and names of atoms in the system
        were different, we could still identify each atom with the
        corresponding atom in the template system.

    The tables that implement this are:

    system
        The list of all the systems in the database

    configuration
        The list of all configurations of all systems, labeled by the
        system they belong to.

    atom
        The list of atoms in all systems and configurations. The
        attributes of the atoms do not depend on the configuration,
        i.e. are unchanging.

    bond
        The list of all bonds in all systems and configurations,
        giving the two atoms that are bonded plus the bond order.

    coordinates
        The fractional or Cartesian coordinates of the atoms as well
        as any other atomic properties that vary by configuration.

    subset
        The instances of subsets.

    element
        The periodic table, used as a foreign key to identify atoms.

    template
        The templates -- labels -- for subsets, which may connect the
        subset to a template system.

    configuration_subset
        A joining table used to define which subsets are "in" a
        configuration.

    subset_atom
        A joining table to define with atoms are in a subset. There is
        an optional field for the template atom if the template has an
        associated system. In this case, the template atom identifies
        which atom in the template system is the same as the given
        atom.

    symmetry
        Information about point or space group symmetry.

    symmetryoperation
        The detailed symmetry operations for the point or space
        groups.

    cell
        The information about the periodicity of a configuration, if
        needed.

    atomset
        A set of atoms, used with a joining table to connect atoms
        with configurations.

    atomset_atom
        The joining table connecting atomsets with their atoms.

    bondset
        A set of bonds, used with a joining table to connect bonds
        with configurations.

    bondset_bond
        The joining table connecting bondsets with their bonds.

    :meta public:
    """

    def __init__(self, parent=None, logger=logger, **kwargs):
        self._parent = parent
        self.logger = logger
        self._attached = {}
        self._current_system_id = None  # The id of the current system
        self._checkpoints = []
        self._filename = None
        self._db = None
        self._cursor = None
        self._items = {}

        if "filename" in kwargs:
            self.filename = kwargs.pop("filename")
        else:
            self.filename = "seamm.db"

    def __del__(self):
        """Destructor: need to close the database if any."""

        # Delete any cached objects
        del self._items

        # And close the database
        if self._db is not None:
            self.db.commit()
            self.db.close()

    def __enter__(self):
        self.db.commit()

        backup = self.parent.copy_system(self, temporary=True)
        self._checkpoints.append(backup)
        return self

    def __exit__(self, etype, value, traceback):
        backup = self._checkpoints.pop()

        if etype is None:
            self.db.commit()

            # Log the changes
            diffs = self.diff(backup)
            if len(diffs) > 0:
                self.version = self.version + 1

            # Not sure why this commit is needed...
            self.db.commit()
        else:
            self.parent.overwrite(self, backup)

        # and delete the copy
        del self.parent[backup.nickname]

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
        table : Table
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

    def __contains__(self, table):
        """Return a boolean indicating if a key exists."""
        # Normal the tablename is used as an identifier, so is quoted with ".
        # Here we need it as a string literal so strip any quotes from it.
        if "." in table:
            schema, table = table.split(".")
            schema = schema.strip('"')
        else:
            schema = "main"

        table = table.strip('"')
        self.cursor.execute(
            "SELECT COUNT(*)"
            f" FROM {schema}.sqlite_master"
            f" WHERE type = 'table' and name = '{table}'"
        )
        return self.cursor.fetchone()[0] == 1

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        # LGTM doesn't like using 'is', but this is equivalent.
        if self.id() == other.id():
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
        if name != other_name and not self.is_attached(other):
            # Attach the other system in order to do comparisons.
            self.attach(other)
            detach = True

        # Check the tables in both systems
        result = True
        for table in in_common:
            if self[table] != other[table]:
                result = False
                break

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    @property
    def cursor(self):
        """A database cursor."""
        return self._cursor

    @property
    def db(self):
        """The database connection."""
        return self._db

    @property
    def db_version(self):
        """The version string for the database."""
        self.cursor.execute("SELECT value FROM metadata WHERE key = 'version'")
        return self.cursor.fetchone()[0]

    @property
    def filename(self):
        """The name of the file (or URI) for the database."""
        return self._filename

    @filename.setter
    def filename(self, value):
        if value != self._filename:
            if self._db is not None:
                self.cursor.close()
                self._db.commit()
                self._db.close()
                self._db = None
                self._cursor = None
        self._filename = value
        if self._filename is not None:
            if self._filename[0:5] == "file:":
                self._db = sqlite3.connect(self._filename, uri=True)
            else:
                self._db = sqlite3.connect(self._filename)
            self._db.row_factory = sqlite3.Row
            self._db.execute("PRAGMA foreign_keys = ON")
            self._cursor = self._db.cursor()
            self._initialize()

    @property
    def names(self):
        """The names of the system."""
        return [row[0] for row in self.db.execute("SELECT name FROM system")]

    @property
    def n_systems(self):
        """The number of systems in the database."""
        return self["system"].n_rows

    @property
    def n_templates(self):
        """The number of templates."""
        return self.templates.n_rows

    @property
    def parent(self):
        """The parent of this, i.e. a Systems object."""
        return self._parent

    @property
    def properties(self):
        """The class to handle the properties."""
        return _Properties(self)

    @property
    def system(self):
        """The current system object."""
        if self.n_systems == 0:
            return None
        if self._current_system_id is None:
            self._current_system_id = self.system_ids[-1]
        return _System(self, self._current_system_id)

    @system.setter
    def system(self, value):
        # Might be a system or an integer
        if isinstance(value, _System):
            value = value.id

        if value not in self.system_ids:
            raise KeyError(f"system '{value}' does not exist.")
        self._current_system_id = value

    @property
    def systems(self):
        """The list of system objects."""
        return [_System(self, _id) for _id in self.system_ids]

    @property
    def system_ids(self):
        """The list of system ids."""
        return [row[0] for row in self.db.execute("SELECT id FROM system")]

    @property
    def templates(self):
        """The defined templates."""
        return _Templates(self)

    def attach(self, other):
        """Attach another system to this one's database.

        Parameters
        ----------
        other : SystemDB
            The other SystemDB object containing the database

        Returns
        -------
        name : str
            The attachment name.
        """
        if not self.is_attached(other):
            # Get a unique name for attaching as
            n = len(self._attached)
            attached_name = f"db_{n}"
            while attached_name in self._attached.values():
                n += 1
                attached_name = f"db_{n}"

            # and attach
            self.db.execute(f"ATTACH DATABASE '{other.filename}' AS '{attached_name}'")
            self._attached[other.filename] = attached_name
        return self._attached[other.filename]

    def attached_as(self, other):
        """The attachment name for another system.

        Parameters
        ----------
        other : SystemDB
            The other SystemDB object containing the attached database

        Returns
        -------
        name : str
            The attachment name.
        """
        return self._attached[other.filename]

    def attributes(self, tablename: str):
        """The attributes -- columns -- of a given table.

        Parameters
        ----------
        tablename : str
            The name of the table, optionally including the schema followed by
            a dot.

        Returns
        -------
        attributes : Dict[str, Any]
            A dictionary of dictionaries for the attributes and their
            descriptors
        """
        if "." in tablename:
            schema, tablename = tablename.split(".")
            sql = f"PRAGMA {schema}.table_info('{tablename}')"
        else:
            sql = f"PRAGMA table_info('{tablename}')"

        result = {}
        for line in self.db.execute(sql):
            result[line["name"]] = {
                "type": line["type"],
                "notnull": bool(line["notnull"]),
                "default": line["dflt_value"],
                "primary key": bool(line["pk"]),
            }
        return result

    def close(self):
        """Close the database."""
        self.filename = None

    def create_system(self, name="", make_current=True):
        """Add a new system.

        Parameters
        ----------
        name : str = None
            A user-friendly name for the system, defaults to no name.
        make_current : bool = True
            If True, make this the current system.

        Returns
        -------
        _System
            The newly created system.
        """

        _id = self["system"].append(name=name)[0]

        if make_current:
            self._current_system_id = _id

        return _System(self, _id)

    def create_table(self, name, cls=_Table, other=None):
        """Create a new table with the given name.

        Parameters
        ----------
        name : str
            The name of the new table.

        cls : Table subclass
            The class of the new table, defaults to Table

        Returns
        -------
        table : class Table
            The new table
        """
        if name in self:
            raise KeyError(f"'{name}' already exists in the system.")

        self._items[name] = cls(self, name, other)
        return self._items[name]

    def delete_system(self, system):
        """Delete an existing system.

        Parameters
        ----------
        system : int or _System
            The system to delete.

        Returns
        -------
        None
        """

        if isinstance(system, _System):
            sid = system.id
        else:
            sid = system

        sql = "DELETE FROM system WHERE id = ?"
        self.db.execute(sql, (sid,))

        if self._current_system_id == sid:
            self._current_system_id = None

    def detach(self, other):
        """Detach an attached system.

        Parameters
        ----------
        other : SystemDB
            The other SystemDB object containing the database
        """
        if self.is_attached(other):
            attached_name = self.attached_as(other)
            self.cursor.execute(f'DETACH DATABASE "{attached_name}"')
            del self._attached[other.filename]

    def diff(self, other):
        """Differences between this system and another."""
        result = {}

        if self.id() == other.id():
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
        name = self.filename
        other_name = other.filename
        detach = False
        if name != other_name and not self.is_attached(other):
            # Attach the other system in order to do comparisons.
            self.attach(other)
            detach = True

        # Check the tables in both systems
        for table in in_common:
            tmp = self[table].diff(other[table])
            if len(tmp) > 0:
                result[f"table '{table}' diffs"] = tmp

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    def find_configurations(self, atomset=None, bondset=None):
        """Return the configurations that have given atom- or bondsets

        Parameters
        ----------
        atomset : int = None
            The id of the atomset.
        bondset : int = None
            The id of the bondset.

        Returns
        -------
        [_Configuration]
        """

        if atomset is not None:
            if bondset is not None:
                self.cursor.execute(
                    "SELECT id FROM configuration WHERE atomset = ? AND bondset = ?",
                    (atomset, bondset),
                )
            else:
                self.cursor.execute(
                    "SELECT id FROM configuration WHERE atomset = ?", (atomset,)
                )
            return [
                _Configuration(_id=cid, system_db=self)
                for cid in self.cursor.fetchall()
            ]
        elif bondset is not None:
            self.cursor.execute(
                "SELECT id FROM configuration WHERE bondset = ?", (bondset,)
            )
            return [
                _Configuration(_id=cid, system_db=self)
                for cid in self.cursor.fetchall()
            ]
        else:
            raise RuntimeError("Must give atomset or bondset.")

    def get_configuration(self, cid):
        """Return the specified configuration.

        Parameters
        ----------
        cid : int
            The id of the configuration.

        Returns
        -------
        _Configuration
            The requested configuration.
        """
        return _Configuration(_id=cid, system_db=self)

    def get_system(self, id_or_name):
        """Get the specified system object.

        Parameters
        ----------
        id_or_name : int or str
            The id (int) or name (str) of the system

        Returns
        -------
        _System
            The requested system.

        Raises
        ------
        ValueError
            If the system does not exist, or more than one have the requested
            name.
        """
        if isinstance(id_or_name, str):
            # The name of the system ... find it.
            sql = "SELECT id FROM SYSTEM WHERE name = ?"
            systems = [x[0] for x in self.db.execute(sql, (id_or_name,))]
            if len(systems) == 0:
                raise ValueError(f"The system '{id_or_name}' does not exist.")
            elif len(systems) > 1:
                raise ValueError(f"There is more than one system named '{id_or_name}'")
            else:
                id_or_name = systems[0]
        return _System(self, id_or_name)

    def is_attached(self, other):
        """Return whether another system is attached to this one.

        Parameters
        ----------
        other : SystemDB
            The other SystemDB object containing the database

        Returns
        -------
        bool
            Whether the database is attached.
        """
        return other.filename in self._attached

    def list(self):
        """Return a list of all the tables in the system."""
        result = []
        for row in self.db.execute(
            "SELECT name" "  FROM sqlite_master" " WHERE type = 'table'"
        ):
            result.append(row["name"])
        return result

    def _initialize(self):
        """Initialize the SQLite database.

        The order is a bit tricky, since many tables reference other tables,
        and hence need to be created after the ones that they reference.

        Start with the standalone tables:
        """

        # If the database is initialized, the metadata table exists.
        # In the future we might need to check the version and upgrade
        # older versions, but now at version 1.0 we are all done!
        if "metadata" not in self:
            # metadata, where we store, get the database version
            table = self["metadata"]
            table.add_attribute("key", coltype="str", pk=True)
            table.add_attribute("value", coltype="str")

            table.append(key="version", value="1.0")
            self.db.commit()

            # The element table
            table = self["element"]
            table.add_attribute("atno", coltype="int", pk=True)
            table.add_attribute("symbol", coltype="str", index="unique")
            table.add_attribute("mass", coltype="float")

            # and fill it from the element data
            for symbol, data in molsystem.elements.data.items():
                table.append(
                    symbol=symbol,
                    atno=data["atomic number"],
                    mass=data["atomic weight"],
                )
            self.db.commit()

            # Symmetry information
            table = self["symmetry"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("group", coltype="str")

            table = self["symmetryoperation"]
            table.add_attribute("symmetry", coltype="int", references="symmetry")
            table.add_attribute("symop", coltype="str")

            # Periodic cell information
            table = self["cell"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("a", coltype="float", default=10.0)
            table.add_attribute("b", coltype="float", default=10.0)
            table.add_attribute("c", coltype="float", default=10.0)
            table.add_attribute("alpha", coltype="float", default=90.0)
            table.add_attribute("beta", coltype="float", default=90.0)
            table.add_attribute("gamma", coltype="float", default=90.0)

            # The systems themselves
            table = self["system"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("name", coltype="str", notnull=True, default="default")

            # The atoms, and the sets of atoms
            table = self["atom"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("atno", coltype="int", references="element")

            table = self["atomset"]
            table.add_attribute("id", coltype="int", pk=True)

            table = self["atomset_atom"]
            table.add_attribute("atomset", coltype="int", references="atomset")
            table.add_attribute("atom", coltype="int", references="atom")
            self.db.execute(
                "CREATE INDEX 'idx_atomset_atom_atomset_atom' "
                'ON atomset_atom ("atomset", "atom")'
            )

            # The bonds, and sets of bonds
            table = self["bond"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("i", coltype="int", references="atom")
            table.add_attribute("j", coltype="int", references="atom")
            table.add_attribute("bondorder", coltype="int", default=1)
            # Only need for infinite networked systems. Think about later!
            # table.add_attribute("offset1", coltype="int", default=0)
            # table.add_attribute("offset2", coltype="int", default=0)
            # table.add_attribute("offset3", coltype="int", default=0)

            table = self["bondset"]
            table.add_attribute("id", coltype="int", pk=True)

            table = self["bondset_bond"]
            table.add_attribute("bondset", coltype="int", references="bondset")
            table.add_attribute("bond", coltype="int", references="bond")
            self.db.execute(
                "CREATE INDEX 'idx_bondset_bond_bondset_bond' "
                'ON bondset_bond ("bondset", "bond")'
            )

            # Now we can set up the configurations
            table = self["configuration"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("name", coltype="str")
            table.add_attribute("system", coltype="int", references="system")
            table.add_attribute("version", coltype="int", notnull=True, default=0)
            table.add_attribute("periodicity", coltype="int", default=0)
            table.add_attribute("coordinate_system", coltype="str", default="Cartesian")
            table.add_attribute("symmetry", coltype="int", references="symmetry")
            table.add_attribute("cell", coltype="int", references="cell")
            table.add_attribute("atomset", coltype="int", references="atomset")
            table.add_attribute("bondset", coltype="int", references="bondset")

            # And coordinates, which depend on configurations
            table = self["coordinates"]
            table.add_attribute(
                "configuration", coltype="int", references="configuration"
            )
            table.add_attribute("atom", coltype="int", references="atom")
            table.add_attribute("x", coltype="float")
            table.add_attribute("y", coltype="float")
            table.add_attribute("z", coltype="float")
            self.db.execute(
                "CREATE INDEX 'idx_coordinates_atom' ON coordinates (\"atom\")"
            )
            self.db.execute(
                "CREATE INDEX idx_coordinates_atom_configuration "
                "    ON coordinates(atom, configuration)"
            )
            # The definition of the subsets -- templates
            table = self["template"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute("name", coltype="str")
            table.add_attribute("canonical_smiles", coltype="str")
            table.add_attribute("category", coltype="str", default="general")
            table.add_attribute(
                "configuration", coltype="int", references="configuration"
            )
            self.db.execute(
                "CREATE UNIQUE INDEX 'idx_template_name_type'"
                '    ON template ("name", "type")'
            )

            # The subsets
            table = self["subset"]
            table.add_attribute("id", coltype="int", pk=True)
            table.add_attribute(
                "configuration", coltype="int", references="configuration"
            )
            table.add_attribute("template", coltype="int", references="template")
            self.db.execute(
                "CREATE INDEX subset_idx_template_configuration "
                "ON subset(template, configuration)"
            )

            # The connection between subsets and the atoms in the system
            table = self["subset_atom"]
            table.add_attribute("atom", coltype="int", references="atom")
            table.add_attribute("subset", coltype="int", references="subset")
            table.add_attribute("templateatom", coltype="int", references="atom")
            self.db.execute(
                "CREATE INDEX 'idx_subset_atom_subset_atom' "
                'ON subset_atom ("subset", "atom")'
            )

            #####################################################
            # The star schema for storing properties and features
            #####################################################

            properties = self.properties
            properties.create_schema()

            self.db.commit()

    def system_exists(self, id_or_name):
        """See if the given system exists.

        Parameters
        ----------
        id_or_name : int or str
            The id (int) or name (str) of the system

        Returns
        -------
        bool
            Whether it exists.
        """
        if isinstance(id_or_name, str):
            # The name of the system ... find it.
            sql = "SELECT COUNT(*) FROM SYSTEM WHERE name = ?"
            self.cursor.execute(sql, (id_or_name,))
        else:
            sql = "SELECT COUNT(*) FROM SYSTEM WHERE id = ?"
            self.cursor.execute(sql, (id_or_name,))

        return self.cursor.fetchone()[0] != 0

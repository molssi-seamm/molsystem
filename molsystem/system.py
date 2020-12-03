# -*- coding: utf-8 -*-

"""A dictionary-like object for holding a system
"""

import collections.abc
from collections import Counter
from functools import reduce
import logging
import math
import sqlite3
from typing import Any, Dict

from molsystem.elemental_data import element_data
from molsystem.table import _Table as Table
from molsystem.atoms import _Atoms as Atoms
from molsystem.subset import _Subsets as Subsets
from molsystem.template import _Template as Template
from molsystem.templateatoms import _Templateatoms as Templateatoms
from molsystem.templatebonds import _Templatebonds as Templatebonds
from molsystem.bonds import _Bonds as Bonds
from molsystem.cell_parameters import _CellParameters as CellParameters

from molsystem.cif import CIFMixin
from molsystem.molfile import MolFileMixin
from molsystem.pdb import PDBMixin
from molsystem.smiles import SMILESMixin
from molsystem.topology import TopologyMixin

logger = logging.getLogger(__name__)


class _System(
    PDBMixin, MolFileMixin, CIFMixin, SMILESMixin, TopologyMixin,
    collections.abc.MutableMapping
):
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
    templateatom -- holds the atoms in a template, if present.
    templatebond -- holds the bonds between the template atoms, if present.

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
    """

    def __init__(self, parent, nickname=None, **kwargs):
        self._parent = parent
        self._nickname = nickname
        self._attached = []
        self._id = 1  # The id of the system. Currently there is only 1
        self._current_configuration = None  # The current configuration
        self._configurations = {}  # template and subset all for configs
        self._checkpoints = []
        self._filename = None
        self._db = None
        self._cursor = None
        self._items = {}
        self._symbol_to_atno = {}
        self._atno_to_symbol = {}
        self._symbol_to_mass = {}
        self._atno_to_mass = {}

        if 'filename' in kwargs:
            self.filename = kwargs.pop('filename')
        else:
            self.filename = 'seamm.db'

    def __del__(self):
        """Destructor: need to close the database if any."""
        if self._db is not None:
            self.db.commit()
            self.db.close()

    def __enter__(self) -> Any:
        self.db.commit()

        backup = self.parent.copy_system(self, temporary=True)
        self._checkpoints.append(backup)
        return self

    def __exit__(self, etype, value, traceback) -> None:
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
        the self._items array to sotre the instances.

        Parameters
        ----------
        key : str
            The table name or name of a multi-table item like atoms.

        Returns
        -------
        table : Table or Atom, Template, Templateatoms, ...
        """
        if key == 'atom' or key == 'atoms':
            if 'atom' not in self._items:
                self._items['atom'] = Atoms(self)
            return self._items['atom']
        elif key == 'template' or key == 'templates':
            if 'template' not in self._items:
                self._items['template'] = Template(self)
            return self._items['template']
        elif key == 'templateatom' or key == 'templateatoms':
            if 'templateatom' not in self._items:
                self._items['templateatom'] = Templateatoms(self)
            return self._items['templateatom']
        elif key == 'templatebond' or key == 'templatebonds':
            if 'templatebond' not in self._items:
                self._items['templatebond'] = Templatebonds(self)
            return self._items['templatebond']
        elif key == 'bond' or key == 'bonds':
            if 'bond' not in self._items:
                self._items['bond'] = Bonds(self)
            return self._items['bond']
        elif key == 'subset' or key == 'subsets':
            if 'subset' not in self._items:
                self._items['subset'] = Subsets(self)
            return self._items['subset']
        elif key == 'cell':
            if 'cell' not in self._items:
                self._items['cell'] = CellParameters(self)
            return self._items['cell']
        else:
            if key not in self._items:
                self._items[key] = Table(self, key)
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
            "SELECT COUNT(*)"
            "  FROM sqlite_master"
            " WHERE type = 'table'"
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
        removed = other_tables - tables
        in_common = tables & other_tables

        if len(added) > 0:
            return False
        if len(removed) > 0:
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
            if Table(self, table) != Table(other, table):
                result = False
                break

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    @property
    def atoms(self):
        """The atoms, which are held as a dictionary of arrays"""
        return self['atom']

    @property
    def bonds(self):
        """The bonds, which are held as a dictionary of arrays"""
        return self['bond']

    @property
    def cell(self):
        """The periodic cell."""
        return self['cell']

    @property
    def configurations(self):
        """The dictionary of configurations."""
        return self._configurations

    @property
    def coordinate_system(self):
        """The coordinates system used, 'fractional' or 'Cartesian'"""
        self.cursor.execute(
            "SELECT coordinatesystem FROM system WHERE id = ?", (self._id,)
        )
        return self.cursor.fetchone()[0]

    @coordinate_system.setter
    def coordinate_system(self, value):
        if value.lower()[0] == 'f':
            self.cursor.execute(
                "UPDATE system SET coordinatesystem = 'fractional'"
                " WHERE id = ?", (self._id,)
            )
        else:
            self.cursor.execute(
                "UPDATE system SET coordinatesystem = 'Cartesian'"
                " WHERE id = ?", (self._id,)
            )

    @property
    def current_configuration(self):
        """The current configuration to work with."""
        return self._current_configuration

    @current_configuration.setter
    def current_configuration(self, value):
        if value not in self.configurations:
            raise KeyError(f"configuration '{value}' doe not exist.")
        self._current_configuration = value

    @property
    def cursor(self):
        return self._cursor

    @property
    def db(self):
        return self._db

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
                self._db = sqlite3.connect(self._filename)
                self._db.row_factory = sqlite3.Row
                self._db.execute('PRAGMA foreign_keys = ON')
                self._cursor = self._db.cursor()
                self._initialize()

    @property
    def name(self):
        """Return the name of this system."""
        self.cursor.execute(
            "SELECT name FROM system WHERE id = ?", (self._id,)
        )
        result = self.cursor.fetchone()
        if result is None:
            return None
        else:
            return result[0]

    @name.setter
    def name(self, value):
        self.cursor.execute(
            "UPDATE system SET name = ? WHERE id = ?", (value, self._id)
        )

    @property
    def nickname(self):
        """The name used in Systems for this system"""
        return self._nickname

    @property
    def n_configurations(self):
        """The number of configurations of the system."""
        return self['configuration'].n_rows

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        self.cursor.execute(
            "SELECT periodicity FROM system WHERE id=?", (self._id,)
        )
        return self.cursor.fetchone()[0]

    @property
    def parent(self):
        """The parent of this, i.e. a Systems object."""
        return self._parent

    @periodicity.setter
    def periodicity(self, value):
        if value < 0 or value > 3:
            raise ValueError('The periodicity must be between 0 and 3.')
        self.cursor.execute(
            "UPDATE system SET periodicity = ? WHERE id = ?",
            (value, self._id)
        )

    @property
    def subsets(self):
        """The subsets"""
        return self['subset']

    @property
    def templates(self):
        """The templates"""
        return self['template']

    @property
    def templateatoms(self):
        """The template atoms"""
        return self['templateatom']

    @property
    def templatebonds(self):
        """The template bonds"""
        return self['templatebond']

    @property
    def version(self):
        """The version of the system, incrementing from 0"""
        self.cursor.execute("SELECT version FROM system")
        return int(self.cursor.fetchone()[0])

    @version.setter
    def version(self, value):
        self.cursor.execute("UPDATE system SET version = ?", (str(value),))
        self.db.commit()

    def add_configuration(
        self,
        system=1,
        name=None,
        symmetry=None,
        cell=None,
        changed_atoms=False,
        changed_bonds=False
    ):
        """Add a new configuration to the system.

        Parameters
        ----------
        system : int = 1
            The system for which this is a configuration. (optional)
        name : str = None
            A textual name for the configuration (optional)
        symmetry : int or str = None
            The id or name of the point or space group (optional)
        cell : Cell or 6-vector = None
            The cell parameters, default to last ones.
        changed_atoms : bool = False
            Whether the atoms have changed in number or identity from the
            previous configuration.
        changed_bonds : bool = False
            Whether the bonding has changed from the previous configuration.

        Returns
        -------
        configuration : int
            The id of the new configuration.

        What needs to be done depends on whether the atoms or bonds are
        changing:

        Case 1: The bonds are changing.
            * create a new template
            * create a subset for the new template
            * if the atoms aren't changing copy the previous connections in
              subset_atom to this instance.
            This routine does not create any templateatoms or bonds.

        Case 2: The atoms are changing
            * create new subset using the previous template
            This routine does not connect the atoms to the new 'all' subset,
            nor does it add any templateatoms or bonds to the new template.

        Case 3: Neither the atoms nor bonding is changing
            * Connect the previous 'all' subset to the configuration in the
              configuration_subset table.
        """
        # Work out the subset and template for 'all'
        if len(self._configurations) == 0:
            changed_bonds = True
            changed_atoms = True
            tid = self['template'].append(name='all', type='all')[0]
            sid = self['subset'].create(template=tid)
        else:
            last_configuration = max(self._configurations)
            last_sid, last_tid = self._configurations[last_configuration]
            sid = last_sid
            tid = last_tid

            if changed_bonds:
                # Case 1
                # If the bonding changed, need a new template and new subset
                tid = self['template'].append(name='all', type='all')
                sid = self['subset'].append(template=tid)
                if not changed_atoms:
                    atom_ids = self['atom'].atoms(
                        configuration=last_configuration
                    )
                    self['subset_atom'].append(subset=sid, atom=atom_ids)
            elif changed_atoms:
                # Case 2
                sid = self['subset'].append(template=tid)[0]
            else:
                # Case 3
                pass

        cid = self['configuration'].append(
            system=system, name=name, symmetry=symmetry
        )[0]

        self['configuration_subset'].append(configuration=cid, subset=sid)
        self._configurations[cid] = (sid, tid)

        return cid

    def all_subset(self, configuration=None):
        if configuration is None:
            configuration = self.current_configuration
        return self._configurations[configuration][0]

    def all_template(self, configuration=None):
        if configuration is None:
            configuration = self.current_configuration
        return self._configurations[configuration][1]

    def attach(self, other):
        """Attach another system to this one's database."""
        if self.is_attached(other.nickname):
            return
        self.db.execute(
            f"ATTACH DATABASE '{other.filename}' AS '{other.nickname}'"
        )
        self._attached.append(other.nickname)

    def is_attached(self, name):
        """Return whether another system is attached to this one."""
        return name in self._attached

    def detach(self, other):
        """Detach an attached system."""
        if self.is_attached(other.name):
            self.cursor.execute(f'DETACH DATABASE "{other.name}"')
            self._attached.remove(other.name)

    def clear(self, configuration=None) -> int:
        """Remove everything from the configuration

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration. Not used if the subset is given.

        Returns
        -------
        int
            Number of atoms
        """
        # Delete the atoms
        self.atoms.remove(configuration=configuration)

        # Delete the template atoms.
        self.templateatoms.remove(template=self.all_template(configuration))

    def create_table(self, name, cls=Table, other=None):
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

    def diff(self, other):
        """Differences between this system and another."""
        result = {}

        if self is other:
            return result

        tables = set(self.list())
        other_tables = set(other.list())

        added = tables - other_tables
        removed = other_tables - tables
        in_common = tables & other_tables

        if len(added) > 0:
            result['tables added'] = list(added)
        if len(removed) > 0:
            result['tables removed'] = list(removed)

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
            tmp = Table(self, table).diff(Table(other, table))
            if len(tmp) > 0:
                result[f"table '{table}' diffs"] = tmp

        # Detach the other database if needed
        if detach:
            self.detach(other)

        return result

    def formula(self, configuration=None):
        """Return the chemical formula of the configuration.

        Returns a tuple with the formula, empirical formula and number of
        formula units (Z).

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.

        Returns
        -------
        formulas : (str, str, int)
            The chemical formula, empirical formula and Z.
        """
        counts = Counter(self.atoms.symbols(configuration))

        # Order the elements ... Merck CH then alphabetical,
        # or if no C or H, then just alphabetically
        formula_list = []
        if 'C' in counts and 'H' in counts:
            formula_list.append(('C', counts.pop('C')))
            formula_list.append(('H', counts.pop('H')))

        for element in sorted(counts.keys()):
            formula_list.append((element, counts[element]))

        counts = []
        for _, count in formula_list:
            counts.append(count)

        formula = []
        for element, count in formula_list:
            if count > 1:
                formula.append(f'{element}{count}')
            else:
                formula.append(element)

        # And the empirical formula
        Z = reduce(math.gcd, counts)
        empirical_formula_list = []
        for element, count in formula_list:
            empirical_formula_list.append((element, int(count / Z)))

        empirical_formula = []
        for element, count in empirical_formula_list:
            if count > 1:
                empirical_formula.append(f'{element}{count}')
            else:
                empirical_formula.append(element)

        return formula, empirical_formula, Z

    def list(self):
        """Return a list of all the tables in the system."""
        result = []
        for row in self.db.execute(
            "SELECT name"
            "  FROM sqlite_master"
            " WHERE type = 'table'"
        ):
            result.append(row['name'])
        return result

    def n_atoms(self, subset=None, configuration=None) -> int:
        """The number of atoms in a subset or configuration

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
        int
            Number of atoms
        """
        return self.atoms.n_atoms(subset=subset, configuration=configuration)

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
        return self.bonds.n_bonds(subset=subset, configuration=configuration)

    def to_atnos(self, symbols):
        """Convert element symbols to atomic numbers.

        Parameters
        ----------
        symbols : [str]
            The atomic symbols

        Returns
        -------
        atnos : [int]
            The corresponding atomic numbers (1..118)
        """
        return [self._symbol_to_atno[x] for x in symbols]

    def to_symbols(self, atnos):
        """Convert atomic numbers to element symbols.

        Parameters
        ----------
        atnos : [int]
            The atomic numbers (1..118)

        Returns
        -------
        symbols : [str]
            The corresponding atomic symbols
        """
        return [self._atno_to_symbol[x] for x in atnos]

    def default_masses(self, symbols=None, atnos=None):
        """Get the atomic mass given atomic symbols or numbers.

        Parameters
        ----------
        symbols : [str] = None
            The atomic symbols
        atnos : [int] = None
            The atomic numbers (1..118)

        Returns
        -------
        masses : [float]
            The default atomic masses
        """
        if symbols is not None:
            return [self._symbol_to_mass[x] for x in symbols]
        if atnos is not None:
            return [self._atno_to_mass[x] for x in atnos]
        else:
            # return all the masses, in order
            return [self._atno_to_mass[x] for x in range(1, 118)]

    def mass(self, subset=None, configuration=None):
        """Return the total atomic masses for the subset or configuration

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
        float
            The summed atomic masses.
        """
        masses = self.atoms.atomic_masses(
            subset=subset, configuration=configuration
        )
        return sum(masses)

    def volume(self, configuration=None):
        """Return the volume of a configuration

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        float
            The volume of the cell.
        """
        if self.periodicity != 3:
            raise RuntimeError('Density is only defined for 3-D systems.')

        return self.cell.cell(configuration=configuration).volume

    def density(self, configuration=None):
        """Return the density of the system.

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        float
            The density of the cell.
        """
        if self.periodicity != 3:
            raise RuntimeError('Density is only defined for 3-D systems.')
        if configuration is None:
            configuration = self.current_configuration
        volume = self.volume(configuration=configuration)
        mass = self.mass(configuration=configuration)

        # converting from g/mol / Å^3 to g/cm^3
        return (mass / volume) * (1.0e+24 / 6.02214076E+23)

    def _initialize(self):
        """Initialize the SQLite database."""
        if 'element' not in self:
            self._initialize_elements()

        if 'symmetry' not in self:
            self._initialize_symmetry()

        if 'cell' not in self:
            self._initialize_cell()

        if 'system' not in self:
            self._initialize_system()
            self.name = self._nickname
        self._id = 1

        if 'configuration' not in self:
            self._initialize_configurations()
        else:
            # Get all the configurations from the database
            for row in self.db.execute(
                "SELECT configuration, subset, template FROM "
                "       configuration_subset, subset, template"
                " WHERE template.name = 'all' AND template.type = 'all'"
                "   AND template.id = template AND subset.id = subset"
            ):
                config = row['configuration']
                self._configurations[config] = (row['subset'], row['template'])
        if 'atom' not in self:
            self._initialize_atoms()

        if 'subset' not in self:
            self._initialize_subsets()

        # If needed, set up the first configuration, and the 'all' subset
        if self.n_configurations == 0:
            self.current_configuration = self.add_configuration()

    def _initialize_system(self):
        """Set up the table for the system."""
        table = self['system']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute(
            'name', coltype='str', notnull=True, default='default'
        )
        table.add_attribute('version', coltype='int', notnull=True, default=0)
        table.add_attribute(
            'periodicity', coltype='int', notnull=True, default=0
        )
        table.add_attribute(
            'coordinatesystem',
            coltype='str',
            notnull=True,
            default='Cartesian'
        )

        table.append(
            id=1,
            name='default',
            version=0,
            periodicity=0,
            coordinatesystem='Cartesian'
        )
        self.db.commit()

    def _initialize_atoms(self):
        """Set up the tables for atoms."""
        table = Table(self, 'atom')
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('atno', coltype='int', references='element')

        table = self['coordinates']
        table.add_attribute(
            'configuration', coltype='int', references='configuration'
        )
        table.add_attribute('atom', coltype='int', references='atom')
        table.add_attribute('x', coltype='float')
        table.add_attribute('y', coltype='float')
        table.add_attribute('z', coltype='float')

    def _initialize_cell(self):
        """Set up the tables for the cell."""
        table = self['cell']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('a', coltype='float', default=10.0)
        table.add_attribute('b', coltype='float', default=10.0)
        table.add_attribute('c', coltype='float', default=10.0)
        table.add_attribute('alpha', coltype='float', default=90.0)
        table.add_attribute('beta', coltype='float', default=90.0)
        table.add_attribute('gamma', coltype='float', default=90.0)

    def _initialize_elements(self):
        """Set up the table of elements."""
        table = self['element']
        table.add_attribute('atno', coltype='int', pk=True)
        table.add_attribute('symbol', coltype='str', index='unique')
        table.add_attribute('mass', coltype='float')

        for symbol, data in element_data.items():
            table.append(
                symbol=symbol,
                atno=data['atomic number'],
                mass=data['atomic weight']
            )
            self._symbol_to_atno[symbol] = data['atomic number']
            self._atno_to_symbol[data['atomic number']] = symbol
            self._symbol_to_mass[symbol] = data['atomic weight']
            self._atno_to_mass[data['atomic number']] = data['atomic weight']
        self.db.commit()

    def _initialize_configurations(self):
        """Set up the table of configurations."""
        table = self['configuration']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('system', coltype='int', references='system')
        table.add_attribute('name', coltype='str')
        table.add_attribute('symmetry', coltype='int', references='symmetry')
        table.add_attribute('cell', coltype='int', references='cell')
        self.db.commit()

    def _initialize_subsets(self):
        """Set up the tables for handling subsets.
        """
        # The definition of the subsets -- templates
        table = self['template']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('name', coltype='str')
        table.add_attribute('type', coltype='str', default='general')
        self.db.execute(
            "CREATE UNIQUE INDEX 'idx_template_name_type'"
            '    ON template ("name", "type")'
        )

        # The subsets
        table = self['subset']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('template', coltype='int', references='template')

        # The template atoms (if any) for the subset
        table = Table(self, 'templateatom')
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('template', coltype='int', references='template')
        table.add_attribute('name', coltype='str')
        table.add_attribute('atno', coltype='int', references='element')

        # The template coordinates (if any) for the subset
        table = Table(self, 'templatecoordinates')
        table.add_attribute(
            'templateatom', coltype='int', references='templateatom'
        )
        table.add_attribute('x', coltype='float')
        table.add_attribute('y', coltype='float')
        table.add_attribute('z', coltype='float')

        # And bonding for the template
        table = self['templatebond']
        table.add_attribute('i', coltype='int', references='templateatom')
        table.add_attribute('j', coltype='int', references='templateatom')
        table.add_attribute('bondorder', coltype='int', default=1)

        # The connection from the subsets to the configurations
        table = self['configuration_subset']
        table.add_attribute(
            'configuration', coltype='int', references='configuration'
        )
        table.add_attribute('subset', coltype='int', references='subset')

        # The connection between subsets and the atoms in the system
        table = self['subset_atom']
        table.add_attribute('atom', coltype='int', references='atom')
        table.add_attribute('subset', coltype='int', references='subset')
        table.add_attribute(
            'templateatom', coltype='int', references='templateatom'
        )

    def _initialize_symmetry(self):
        """Set up the tables for symmetry."""
        table = self['symmetry']
        table.add_attribute('id', coltype='int', pk=True)
        table.add_attribute('group', coltype='str')

        table = self['symmetryoperation']
        table.add_attribute('symmetry', coltype='int', references='symmetry')
        table.add_attribute('symop', coltype='str')

    def attributes(self, tablename: str) -> Dict[str, Any]:
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
        if '.' in tablename:
            schema, tablename = tablename.split('.')
            sql = f"PRAGMA {schema}.table_info('{tablename}')"
        else:
            sql = f"PRAGMA table_info('{tablename}')"

        result = {}
        for line in self.db.execute(sql):
            result[line['name']] = {
                'type': line['type'],
                'notnull': bool(line['notnull']),
                'default': line['dflt_value'],
                'primary key': bool(line['pk'])
            }
        return result

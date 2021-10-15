# -*- coding: utf-8 -*-

from collections import Counter
from functools import reduce
import logging
import math

from .atoms import _Atoms
from .bonds import _Bonds
from .cell import _Cell
from .cif import CIFMixin
from .configuration_properties import _ConfigurationProperties
from .molfile import MolFileMixin
from .openbabel import OpenBabelMixin
from .rdkit_ import RDKitMixin
from .pdb import PDBMixin
from .smiles import SMILESMixin
from .subsets import _Subsets
from .symmetry import _Symmetry
from .topology import TopologyMixin

logger = logging.getLogger(__name__)


class _Configuration(
    PDBMixin,
    MolFileMixin,
    CIFMixin,
    SMILESMixin,
    TopologyMixin,
    OpenBabelMixin,
    RDKitMixin,
    object,
):
    """A configuration (conformer) of a system.

    :meta public:
    """

    def __init__(self, _id, system_db):

        self._id = _id
        self._system_db = system_db

        self._name = None
        self._system = None
        self._atomset = None
        self._bondset = None
        self._cell_id = None
        self._coordinate_system = None
        self._periodicity = None
        self._symmetry_id = None

    def __enter__(self):
        """Copy the tables to a backup for a 'with' statement."""
        self.system_db["configuration"].__enter__()
        self.atoms.__enter__()
        self.bonds.__enter__()
        # Avoid periodicity check
        _Cell(self).__enter__()
        return self

    def __exit__(self, etype, value, traceback):
        """Handle returning from a 'with' statement."""
        # Save the version, beacuse all of the subitems will update it.
        version = self.version

        self.system_db["configuration"].__exit__(etype, value, traceback)
        self.atoms.__exit__(etype, value, traceback)
        self.bonds.__exit__(etype, value, traceback)
        _Cell(self).__exit__(etype, value, traceback)

        if etype is None:
            self.version = version + 1

        return False

    @property
    def atoms(self):
        """The atoms for this configuration."""
        return _Atoms(self)

    @property
    def atomset(self):
        """The id of the atom set for this configuration."""
        if self._atomset is None:
            # Cache the atomset for this configuration
            self.cursor.execute(
                "SELECT atomset FROM configuration WHERE id = ?", (self.id,)
            )
            self._atomset = self.cursor.fetchone()[0]
        if self._atomset is None:
            # No atomset, so create one and update this configuration
            self._atomset = self.system_db["atomset"].append(n=1)[0]
            self.db.execute(
                "UPDATE configuration SET atomset = ? WHERE id = ?",
                (self._atomset, self.id),
            )
            self.db.commit()
        return self._atomset

    @atomset.setter
    def atomset(self, value):
        if self._atomset is None:
            # Cache the atomset for this configuration
            self.cursor.execute(
                "SELECT atomset FROM configuration WHERE id = ?", (self.id,)
            )
            tmp = self.cursor.fetchone()
            self._atomset = tmp[0]
        if self._atomset is None:
            # No atomset, so update this configuration
            self.db.execute(
                "UPDATE configuration SET atomset = ? WHERE id = ?",
                (value, self.id),
            )
            self.db.commit()
            self._atomset = value
        elif value == self._atomset:
            # Do nothing
            pass
        else:
            raise RuntimeError("The atomset is already set!")

    @property
    def bonds(self):
        """The bonds for this configuration."""
        return _Bonds(self)

    @property
    def bondset(self):
        """The id of the bond set for this configuration."""
        if self._bondset is None:
            # Cache the bondset for this configuration
            self.cursor.execute(
                "SELECT bondset FROM configuration WHERE id = ?", (self.id,)
            )
            self._bondset = self.cursor.fetchone()[0]
        if self._bondset is None:
            # No bondset, so create one and update this configuration
            self._bondset = self.system_db["bondset"].append(n=1)[0]
            self.db.execute(
                "UPDATE configuration SET bondset = ? WHERE id = ?",
                (self._bondset, self.id),
            )
            self.db.commit()
        return self._bondset

    @bondset.setter
    def bondset(self, value):
        if self._bondset is None:
            # Cache the bondset for this configuration
            self.cursor.execute(
                "SELECT bondset FROM configuration WHERE id = ?", (self.id,)
            )
            self._bondset = self.cursor.fetchone()[0]
        if self._bondset is None:
            # No bondset, so update this configuration
            self.db.execute(
                "UPDATE configuration SET bondset = ? WHERE id = ?",
                (value, self.id),
            )
            self.db.commit()
            self._bondset = value
        elif value == self._bondset:
            # No change
            pass
        else:
            raise RuntimeError("The bondset is already set!")

    @property
    def cell(self):
        """The periodic (unit) cell for this configuration.

        Raises
        ------
        TypeError
            If the system is not periodic.

        Returns
        -------
        _Cell
            The Cell object for this configuration.
        """
        if self.periodicity == 0:
            raise TypeError("The configuration is not periodic!")
        return _Cell(self)

    @property
    def cell_id(self):
        """The id of the cell in the cell table."""
        if self._cell_id is None:
            self.cursor.execute(
                "SELECT cell FROM configuration WHERE id = ?", (self.id,)
            )
            self._cell_id = self.cursor.fetchone()[0]
            if self._cell_id is None:
                self._cell_id = self.system_db["cell"].append(n=1)[0]
                sql = "UPDATE configuration SET cell = ? WHERE id = ?"
                self.db.execute(sql, (self._cell_id, self.id))
                self.db.commit()
        return self._cell_id

    @property
    def coordinate_system(self):
        """The coordinate system for this configuration."""
        if self._coordinate_system is None:
            self.cursor.execute(
                "SELECT coordinate_system FROM configuration WHERE id = ?", (self.id,)
            )
            self._coordinate_system = self.cursor.fetchone()[0]
        return self._coordinate_system

    @coordinate_system.setter
    def coordinate_system(self, value):
        if value.lower()[0] == "f":
            self.cursor.execute(
                "UPDATE configuration SET coordinate_system = 'fractional'"
                " WHERE id = ?",
                (self.id,),
            )
            self._coordinate_system = "fractional"
        else:
            self.cursor.execute(
                "UPDATE configuration SET coordinate_system = 'Cartesian'"
                " WHERE id = ?",
                (self.id,),
            )
            self._coordinate_system = "Cartesian"
        self.db.commit()

    @property
    def cursor(self):
        """The database connection."""
        return self.system_db.cursor

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def density(self):
        """Return the density of this configuration.

        Returns
        -------
        float
            The density of the cell.
        """
        if self.periodicity != 3:
            raise RuntimeError("Density is only defined for 3-D systems.")

        volume = self.volume
        mass = self.mass

        # converting from g/mol / Å^3 to g/cm^3
        return (mass / volume) * (1.0e24 / 6.02214076e23)

    @property
    def formula(self):
        """Return the chemical formula of the configuration.

        Returns a tuple with the formula, empirical formula and number of
        formula units (Z).

        Returns
        -------
        formulas : (str, str, int)
            The chemical formula, empirical formula and Z.
        """
        counts = Counter(self.atoms.symbols)

        # Order the elements ... Hill order if C, CH then alphabetical,
        # or if no C, then just alphabetically
        formula_list = []
        if "C" in counts:
            formula_list.append(("C", counts.pop("C")))
            if "H" in counts:
                formula_list.append(("H", counts.pop("H")))

        for element in sorted(counts.keys()):
            formula_list.append((element, counts[element]))

        counts = []
        for _, count in formula_list:
            counts.append(count)

        formula = []
        for element, count in formula_list:
            if count > 1:
                formula.append(f"{element}{count}")
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
                empirical_formula.append(f"{element}{count}")
            else:
                empirical_formula.append(element)

        return " ".join(formula), " ".join(empirical_formula), Z

    @property
    def id(self):
        """The id of this configuration."""
        return self._id

    @property
    def mass(self):
        """Return the total atomic mass.

        Returns
        -------
        float
            The summed atomic masses.
        """
        masses = self.atoms.atomic_masses

        return sum(masses)

    @property
    def n_atoms(self) -> int:
        """The number of atoms.

        Returns
        -------
        int
            Number of atoms
        """
        return self.atoms.n_atoms

    @property
    def name(self):
        """The name of the configuration."""
        if self._name is None:
            self.cursor.execute(
                "SELECT name FROM configuration WHERE id = ?", (self.id,)
            )
            self._name = self.cursor.fetchone()[0]
        return self._name

    @name.setter
    def name(self, value):
        self.db.execute(
            "UPDATE configuration SET name = ? WHERE id = ?", (value, self.id)
        )
        self.db.commit()
        self._name = value

    @property
    def n_bonds(self) -> int:
        """The number of bonds.

        Returns
        -------
        int
            Number of bonds
        """
        return self.bonds.n_bonds

    @property
    def periodicity(self):
        """The periodicity of the system, 0, 1, 2 or 3"""
        if self._periodicity is None:
            self.cursor.execute(
                "SELECT periodicity FROM configuration WHERE id = ?", (self.id,)
            )
            self._periodicity = self.cursor.fetchone()[0]
        return self._periodicity

    @periodicity.setter
    def periodicity(self, value):
        if value < 0 or value > 3:
            raise ValueError("The periodicity must be between 0 and 3.")
        self.cursor.execute(
            "UPDATE configuration SET periodicity = ? WHERE id = ?", (value, self.id)
        )
        self.db.commit()
        self._periodicity = value

    @property
    def properties(self):
        """The class to handle the properties for this configuration."""
        return _ConfigurationProperties(self)

    @property
    def subsets(self):
        """The subsets"""
        return _Subsets(self)

    @property
    def symmetry(self):
        """The periodic (unit) symmetry for this configuration."""
        return _Symmetry(self.system_db, self.symmetry_id)

    @property
    def symmetry_id(self):
        """The id of the symmetry in the symmetry table."""
        if self._symmetry_id is None:
            self.cursor.execute(
                "SELECT symmetry FROM configuration WHERE id = ?", (self.id,)
            )
            self._symmetry_id = self.cursor.fetchone()[0]
        return self._symmetry_id

    @property
    def system(self):
        """The system that we belong to."""
        if self._system is None:
            self.cursor.execute(
                "SELECT system FROM configuration WHERE id = ?", (self.id,)
            )
            self._system = self.system_db.get_system(self.cursor.fetchone()[0])
        return self._system

    @property
    def system_db(self):
        """The system that we belong to."""
        return self._system_db

    @property
    def version(self):
        """The version of the system, incrementing from 0"""
        self.cursor.execute(
            "SELECT version FROM configuration WHERE id = ?", (self.id,)
        )
        return int(self.cursor.fetchone()[0])

    @version.setter
    def version(self, value):
        self.cursor.execute(
            "UPDATE configuration SET version = ? WHERE id = ?", (value, self.id)
        )
        self.db.commit()

    @property
    def volume(self):
        """Return the volume of this configuration.

        Returns
        -------
        float
            The volume of the cell.
        """
        if self.periodicity != 3:
            raise RuntimeError("Density is only defined for 3-D systems.")

        return self.cell.volume

    def clear(self) -> int:
        """Delete everything from the configuration."""
        # Delete the atoms
        self.atoms.delete("all")

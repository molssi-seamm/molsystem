# -*- coding: utf-8 -*-

from collections import Counter
from functools import reduce
import logging
import math

import seekpath
import spglib

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
from .table import _Table
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
        self._charge = None
        self._coordinate_system = None
        self._spin_multiplicity = 0
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
    def charge(self):
        """The charge of the system, 0, ±1, ±2, ±3, ..."""
        if self._charge is None:
            self.cursor.execute(
                "SELECT charge FROM configuration WHERE id = ?", (self.id,)
            )
            self._charge = self.cursor.fetchone()[0]
        return self._charge

    @charge.setter
    def charge(self, value):
        self.cursor.execute(
            "UPDATE configuration SET charge = ? WHERE id = ?", (value, self.id)
        )
        self.db.commit()
        self._charge = value

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
    def spin_multiplicity(self):
        """The spin_multiplicity of the system, 0, 1, 2, 3, ..."""
        if self._spin_multiplicity == 0:
            self.cursor.execute(
                "SELECT spin_multiplicity FROM configuration WHERE id = ?", (self.id,)
            )
            multiplicity = self.cursor.fetchone()[0]
            if multiplicity is None or multiplicity == 0:
                n_electrons = sum(self.atoms.atomic_numbers) - self.charge
                if n_electrons % 2 == 0:
                    multiplicity = 1
                else:
                    multiplicity = 2
                self.spin_multiplicity = multiplicity
            else:
                self._spin_multiplicity = multiplicity
        return self._spin_multiplicity

    @spin_multiplicity.setter
    def spin_multiplicity(self, value):
        self.cursor.execute(
            "UPDATE configuration SET spin_multiplicity = ? WHERE id = ?",
            (value, self.id),
        )
        self.db.commit()
        self._spin_multiplicity = value

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
    def n_active_electrons(self):
        """The number of active electrons in the configuration"""
        self.cursor.execute(
            "SELECT n_active_electrons FROM configuration WHERE id = ?", (self.id,)
        )
        return self.cursor.fetchone()[0]

    @n_active_electrons.setter
    def n_active_electrons(self, value):
        self.cursor.execute(
            "UPDATE configuration SET n_active_electrons = ? WHERE id = ?",
            (value, self.id),
        )
        self.db.commit()

    @property
    def n_active_orbitals(self):
        """The number of active orbitals in the configuration"""
        self.cursor.execute(
            "SELECT n_active_orbitals FROM configuration WHERE id = ?", (self.id,)
        )
        return self.cursor.fetchone()[0]

    @n_active_orbitals.setter
    def n_active_orbitals(self, value):
        self.cursor.execute(
            "UPDATE configuration SET n_active_orbitals = ? WHERE id = ?",
            (value, self.id),
        )
        self.db.commit()

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
    def state(self):
        """The electronic state of the configuration. Either a number or e.g. 2T1g
        for the second T1g state."""
        self.cursor.execute("SELECT state FROM configuration WHERE id = ?", (self.id,))
        return self.cursor.fetchone()[0]

    @state.setter
    def state(self, value):
        self.cursor.execute(
            "UPDATE configuration SET state = ? WHERE id = ?",
            (value, self.id),
        )
        self.db.commit()

    @property
    def symmetry(self):
        """The periodic (unit) symmetry for this configuration."""
        return _Symmetry(self)

    @property
    def symmetry_id(self):
        """The id of the symmetry in the symmetry table."""
        if self._symmetry_id is None:
            self.cursor.execute(
                "SELECT symmetry FROM configuration WHERE id = ?", (self.id,)
            )
            self._symmetry_id = self.cursor.fetchone()[0]
            if self._symmetry_id is None:
                table = _Table(self.system_db, "symmetry")
                self._symmetry_id = table.append()[0]
                self.cursor.execute(
                    "UPDATE configuration SET symmetry = ? WHERE id = ?",
                    (self._symmetry_id, self.id),
                )
                self.db.commit()
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

    def primitive_cell(self, symprec=1.0e-05, spg=False):
        """Find the symmetry of periodic systems and transform to conventional cell.

        Parameters
        ----------
        symprec : float
            Distance tolerance in Cartesian coordinates to find crystal symmetry.

            For atomic positions, roughly speaking, two position vectors x and x’ in
            Cartesian coordinates are considered to be the same if |x' - x| <
            symprec. For more details, see the spglib paper, Sec. II-A.

            The angle distortion between basis vectors is converted to a length and
            compared with this distance tolerance. For more details, see the spglib
            paper, Sec. IV-A. It is possible to specify angle tolerance explicitly, see
            angle_tolerance.

        spg : bool = False
            Whether to use ``spglib`` or ``seekpath`` (the default)

        Returns
        -------
        ([[float*3]*3], [[float*3]*natoms], [int*natoms])
            A tuple with the cell vectors, fractional coordinates and atomic numbers.
        """
        lattice_in = self.cell.vectors()
        fractionals_in = self.atoms.get_coordinates(fractionals=True)
        atomic_numbers_in = self.atoms.atomic_numbers
        cell_in = (lattice_in, fractionals_in, atomic_numbers_in)

        if spg:
            lattice, fractionals, atomic_numbers = spglib.standardize_cell(
                cell_in, to_primitive=True, no_idealize=True, symprec=symprec
            )
        else:
            dataset = seekpath.get_path(cell_in, symprec=symprec)
            lattice = dataset["primitive_lattice"].tolist()
            fractionals = dataset["primitive_positions"].tolist()
            atomic_numbers = dataset["primitive_types"].tolist()

        return lattice, fractionals, atomic_numbers

    def symmetrize(self, symprec=1.0e-05, angle_tolerance=None, spg=False):
        """Find the symmetry of periodic systems and transform to conventional cell.

        Parameters
        ----------
        symprec : float
            Distance tolerance in Cartesian coordinates to find crystal symmetry.

            For atomic positions, roughly speaking, two position vectors x and x’ in
            Cartesian coordinates are considered to be the same if |x' - x| <
            symprec. For more details, see the spglib paper, Sec. II-A.

            The angle distortion between basis vectors is converted to a length and
            compared with this distance tolerance. For more details, see the spglib
            paper, Sec. IV-A. It is possible to specify angle tolerance explicitly, see
            angle_tolerance.

        angle_tolerance : float = None
            Tolerance of angle between basis vectors in degrees to be tolerated in the
            symmetry finding. If None, then the ``symprec`` above is used. This is the
            recommended approach!

        spg : bool = False
            Whether to use ``spglib`` or ``seekpath`` (the default)
        """
        if self.periodicity == 3:
            lattice_in = self.cell.vectors()
            fractionals_in = self.atoms.get_coordinates(fractionals=True)
            atomic_numbers_in = self.atoms.atomic_numbers
            cell_in = (lattice_in, fractionals_in, atomic_numbers_in)

            logger.debug("Lattice input to symmetrize")
            for vector in lattice_in:
                logger.debug(
                    f"   {vector[0]:7.3f}   {vector[1]:7.3f}   {vector[2]:7.3f}"
                )
            logger.debug("")
            logger.debug("Fractionals")
            for vector, atno in zip(fractionals_in, atomic_numbers_in):
                logger.debug(
                    f"   {atno:2} {vector[0]:7.3f}   {vector[1]:7.3f}   "
                    f"{vector[2]:7.3f}"
                )

            if spg:
                if angle_tolerance is not None:
                    dataset = spglib.get_symmetry_dataset(
                        cell_in, symprec=symprec, angle_tolerance=angle_tolerance
                    )
                else:
                    dataset = spglib.get_symmetry_dataset(cell_in, symprec=symprec)

                lattice = dataset["std_lattice"].tolist()
                fractionals = dataset["std_positions"].tolist()
                atomic_numbers = dataset["std_types"].tolist()
                space_group = dataset["international"].strip("'")
            else:
                if angle_tolerance is not None:
                    dataset = seekpath.get_path(
                        cell_in, symprec=symprec, angle_tolerance=angle_tolerance
                    )
                else:
                    dataset = seekpath.get_path(cell_in, symprec=symprec)

                lattice = dataset["conv_lattice"].tolist()
                fractionals = dataset["conv_positions"].tolist()
                atomic_numbers = dataset["conv_types"].tolist()
                space_group = dataset["spacegroup_international"].strip("'")

            logger.debug("")
            logger.debug("Lattice after symmetrize")
            for vector in lattice:
                logger.debug(
                    f"   {vector[0]:7.3f}   {vector[1]:7.3f}   {vector[2]:7.3f}"
                )
            logger.debug("")
            logger.debug("Fractionals")
            for vector, atno in zip(fractionals, atomic_numbers):
                logger.debug(
                    f"   {atno:2} {vector[0]:7.3f}   {vector[1]:7.3f}   "
                    f"{vector[2]:7.3f}"
                )

            self.cell.from_vectors(lattice)

            self.atoms.delete("all")

            xs = [f[0] for f in fractionals]
            ys = [f[1] for f in fractionals]
            zs = [f[2] for f in fractionals]

            self.coordinate_system = "fractional"
            self.atoms.append(atno=atomic_numbers, x=xs, y=ys, z=zs)
            self.symmetry.group = space_group

    def update(
        self,
        coordinates,
        fractionals=True,
        atomic_numbers=None,
        lattice=None,
        space_group=None,
        symprec=1.0e-05,
        angle_tolerance=None,
        spg=False,
    ):
        """Update the system, checking symmetry.

        Parameters
        ----------
        coordinates : [[3*float]*natoms]
            The coordinates.

        fractionals : bool = True
            Whether fractional or Cartesian coordinates.

        atomic_numbers : [float*natoms]
            The atomic numbers of the atoms

        lattice : [[float*3]*3] = None
            The cell vectors

        space_group : str = None
            The original space group.

        symprec : float
            Distance tolerance in Cartesian coordinates to find crystal symmetry.

            For atomic positions, roughly speaking, two position vectors x and x’ in
            Cartesian coordinates are considered to be the same if |x' - x| <
            symprec. For more details, see the spglib paper, Sec. II-A.

            The angle distortion between basis vectors is converted to a length and
            compared with this distance tolerance. For more details, see the spglib
            paper, Sec. IV-A. It is possible to specify angle tolerance explicitly, see
            angle_tolerance.

        angle_tolerance : float = None
            Tolerance of angle between basis vectors in degrees to be tolerated in the
            symmetry finding. If None, then the ``symprec`` above is used. This is the
            recommended approach!

        spg : bool = False
            Whether to use ``spglib`` or ``seekpath`` (the default)
        """
        text = ""

        cell_in = (lattice, coordinates, atomic_numbers)

        if spg:
            if angle_tolerance is not None:
                dataset = spglib.get_symmetry_dataset(
                    cell_in, symprec=symprec, angle_tolerance=angle_tolerance
                )
            else:
                dataset = spglib.get_symmetry_dataset(cell_in, symprec=symprec)

            lattice_out = dataset["std_lattice"].tolist()
            fractionals_out = dataset["std_positions"].tolist()
            atomic_numbers_out = dataset["std_types"].tolist()
            space_group_out = dataset["international"].strip("'")
        else:
            if angle_tolerance is not None:
                dataset = seekpath.get_path(
                    cell_in, symprec=symprec, angle_tolerance=angle_tolerance
                )
            else:
                dataset = seekpath.get_path(cell_in, symprec=symprec)

            lattice_out = dataset["conv_lattice"].tolist()
            fractionals_out = dataset["conv_positions"].tolist()
            atomic_numbers_out = dataset["conv_types"].tolist()
            space_group_out = dataset["spacegroup_international"].strip("'")

        logger.debug("Lattice")
        for vector in lattice_out:
            logger.debug(f"   {vector[0]:7.3f}   {vector[1]:7.3f}   {vector[2]:7.3f}")
        logger.debug("")
        logger.debug("Fractionals")
        for vector, atno in zip(fractionals_out, atomic_numbers_out):
            logger.debug(
                f"   {atno:2} {vector[0]:7.3f}   {vector[1]:7.3f}   {vector[2]:7.3f}"
            )
        logger.debug(f"{space_group_out}")

        if space_group != space_group_out:
            text = f"The space group changed from {space_group} to {space_group_out}."
            logger.debug(text)

        self.cell.from_vectors(lattice_out)

        self.atoms.delete("all")

        xs = [f[0] for f in fractionals_out]
        ys = [f[1] for f in fractionals_out]
        zs = [f[2] for f in fractionals_out]

        self.coordinate_system = "fractional"
        self.atoms.append(atno=atomic_numbers_out, x=xs, y=ys, z=zs)
        self.symmetry.group = space_group_out

        return text

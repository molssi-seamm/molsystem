# -*- coding: utf-8 -*-

from fractions import Fraction
import logging
import pprint  # noqa: F401
import re

import numpy as np
import spglib

logger = logging.getLogger(__name__)
# logger.setLevel("INFO")
logger.setLevel("WARNING")


class _Symmetry(object):
    """A class to handle point and space group symmetry."""

    spgno_to_hall = None
    spgname_to_hall = None
    spgname_to_system = None
    hall_to_spgname = None
    hall_to_hall_symbol = None
    hall_to_IT_symbol = None

    def __init__(self, configuration):
        """Initialize from the database.

        Parameters
        ----------
        system_db : SystemDB
            The SystemDB instance that we are working with.
        _id : int
            The id of this particular symmetry.
        """
        self._configuration = configuration
        self._system = self._configuration.system
        self._system_db = self._system.system_db
        self._id = configuration.symmetry_id
        self._operators = None  # 4x4 symmetry operation matrices
        self._symop_products = None  # Square matrix of symops that are products of 2
        self._atom_generators = None  # Sym operations that create the symmetric atoms
        self._symop_to_atom = None  # The symmetry atom resulting from symmetry ops
        self._atom_to_asymmetric_atom = None  # The asymmetric atom for each atom
        self._bond_to_asymmetric_bond = None  # The asymmetric bond for each bond
        self._bond_atoms = None  # The pair of symmetric atoms in each bond
        self._bond_offsets = None  # The triplet of cell offsets for each bond

        super().__init__()

    @staticmethod
    def filter_atoms(atoms, coordinates, group=None, sym_ops=None, operators=None):
        """Check and filter the atoms in a periodic system.

        There is an issue with some CIF files, notably those from the Cambridge
        structural database, where there are atoms in the asymmetric unit that
        are duplicates. It appears that the CIF file is written with entire molecules
        in some cases. So if the molecule is on a symmetry element, the atoms that are
        related by symmetry are explicitly in the file.

        To handle this we need to check for duplicates and remove them. Not a bad check
        to have, anyway.

        Parameters
        ----------
        atoms : str or int
            The atomic symbols or numbers of the atoms.
        coordinates : array_like
            The coordinates of the atoms.
        group : str
            The space group symbol.
        sym_ops : array_like
            The symmetry operations as strings.
        operators : array_like
            The symmetry operators.
        """
        if group is not None:
            # Get the symmetry operators
            if sym_ops is not None or operators is not None:
                raise RuntimeError(
                    "Cannot specify both group and sym_ops or operators."
                )
            operators = _Symmetry.get_operators(group)
        elif sym_ops is not None:
            if operators is not None or group is not None:
                raise RuntimeError(
                    "Cannot specify both sym_ops and group or operators."
                )
            operators = _Symmetry.get_operators(sym_ops)
        elif operators is not None:
            if sym_ops is not None or group is not None:
                raise RuntimeError(
                    "Cannot specify both operators and group or sym_ops."
                )
        else:
            raise RuntimeError("Must specify either group or sym_ops or operators.")

        self = _Symmetry()

        symbols = self.configuration.atoms.asymmetric_symbols
        n_atoms = self.configuration.n_asymmetric_atoms
        logger.info(f"Expanding {n_atoms} asymmetric atoms to full cell.")
        if n_atoms == 0:
            self._atom_generators = None
            self._symop_to_atom = None
        self._atom_generators = []  # Symmetry ops that create the symmetric atoms
        self._symop_to_atom = []  # The symmetry atom resulting from symmetry ops
        if self.configuration.periodicity == 3:
            if len(operators) == 1:
                # P1 is a special case. Nothing to do.
                self._atom_generators = [[0] for i in range(n_atoms)]
                self._symop_to_atom = [[i] for i in range(n_atoms)]
                self._atom_to_asymmetric_atom = [*range(n_atoms)]
            else:
                uvw0 = self.configuration.atoms.get_coordinates(
                    as_array=True, asymmetric=True
                )
                if uvw0.shape[0] != n_atoms:
                    raise RuntimeError(
                        f"Mismatch of number of atoms in symmetry: {uvw0.shape[0]} != "
                        f"{n_atoms}"
                    )
                logger.debug(
                    f"""
{n_atoms=}
{uvw0.shape=}"
Original coordinates
{uvw0}
                    """
                )

                uvw = np.ndarray((n_atoms, 4))
                uvw[:, 0:3] = uvw0[:, :]
                logger.debug(f"\nExpanded coordinates\n{uvw}")
                uvw[:, 3] = 1
                # logger.debug(self.symops)
                logger.debug(f"\n{operators.shape=}")
                # logger.debug(operators)
                # logger.debug(f"{uvw.shape=}")
                # logger.debug("Expanded coordinates")
                # logger.debug(uvw)
                xformed = np.einsum("ijk,lk", operators, uvw)
                logger.debug(f"\n{xformed.shape=}\n{xformed}")

                # For comparison, bring all atoms into the cell [0..1)
                tmp = xformed - np.floor(xformed)

                logger.debug(f"\ntmp = Transformed and floored coordinates\n{tmp}")

                # Keep track of atoms and coordinates found so can check for duplicates
                found = {}

                n_sym_atoms = 0
                n_asymmetric_atoms = 0
                for i in range(n_atoms):
                    values, I1, I2 = np.unique(
                        np.round(tmp[:, :3, i], 4),
                        axis=0,
                        return_index=True,
                        return_inverse=True,
                    )
                    # pprint.pprint(np.round(tmp[:, :3, i], 4).tolist())
                    logger.debug(
                        f"""
{i=}: {symbols[i]}
{np.round(tmp[:, :3, i], 4)}
nvalues
{values}
I1
{I1}
I2
{I2}
                        """
                    )
                    # Check for duplicates
                    if tuple(values[0]) in found:
                        if found[tuple(values[0])] != symbols[i]:
                            raise RuntimeError(
                                "Duplicate atoms with different symbols: "
                                f"{symbols[i]} != {found[tuple(values[0])]} at "
                                f"{values[0]}"
                            )
                        else:
                            # Same element, so ignore
                            logger.debug("\nDuplicate atom, ignoring")
                            continue

                    for value in values:
                        found[tuple(value)] = symbols[i]

                    n_asymmetric_atoms += 1
                    I2 += n_sym_atoms
                    self._atom_generators.append(I1.tolist())
                    self._symop_to_atom.append(I2.tolist())
                    n_sym_atoms += I1.shape[0]
                tmp = []
                for i, generators in enumerate(self._atom_generators):
                    tmp.extend([i] * len(generators))
                self._atom_to_asymmetric_atom = tmp
                logger.info(f"{n_asymmetric_atoms=} {n_sym_atoms=}")

    @staticmethod
    def get_operators(self, value, periodicity=3):
        """Get the symmetry operators as 4x4 matrices.

        Parameters
        ----------
        value : str or [str]
            The symmetry group or space group, or symmetry operators as strings.
        periodicity : int
            The periodicity of the structure. 3 for 3D, 2 for 2D, 1 for 1D, or 0.

        Returns
        -------
        operators : np.ndarray
            The symmetry operators as 4x4 matrices.
        """
        if isinstance(value, str):
            if value == "":
                raise RuntimeError("No group specified.")
            if self.periodicity == 0:
                if value != "C1":
                    raise NotImplementedError("Point groups not implemented yet!")
                W4 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

                W4s = []
                W4s.append(W4)
                operators = np.array(W4s, dtype=float)
            else:
                # Get the 4x4 augmented matrices
                hall = _Symmetry.to_hall(value)
                value = _Symmetry.hall_to_spacegroup_name(hall)
                data = spglib.get_symmetry_from_database(hall)
                W4s = []
                for W, w in zip(
                    data["rotations"].tolist(), data["translations"].tolist()
                ):
                    W4 = []
                    for Wrow, w_i in zip(W, w):
                        W4.append([*Wrow, w_i])
                    W4.append([0, 0, 0, 1])
                    W4s.append(W4)
                operators = np.array(W4s, dtype=float)
        elif isinstance(value, list):
            W4s = []
            operators = np.array(W4s, dtype=float)
        return operators

    @staticmethod
    def hall_to_spacegroup_name(hall):
        """Return the International name including setting for a hall number."""
        if _Symmetry.hall_to_spgname is None:
            _Symmetry.spgname_to_hall = None
            _Symmetry.spacegroup_names_to_hall()
        return _Symmetry.hall_to_spgname[hall]

    @staticmethod
    def spacegroup_names_to_hall():
        """Dictionary of Hall number for spacegroup names."""
        if _Symmetry.spgname_to_hall is None:
            system_name = {
                "international_full": "name_H-M_full",
                "international": "name_H-M_alt",
                "international_short": "name_H-M_short",
                "hall_symbol": "name_Hall",
            }
            # Initialize the symmetry data
            _Symmetry.spgname_to_hall = {}
            _Symmetry.hall_to_spgname = {}
            _Symmetry.hall_to_hall_symbol = {}
            _Symmetry.hall_to_IT_number = {}
            if _Symmetry.spgname_to_system is None:
                _Symmetry.spgname_to_system = {}
            for hall in range(1, 530):
                data = spglib.get_spacegroup_type(hall)
                # pprint.pprint(data)
                choice = data["choice"]
                _Symmetry.hall_to_hall_symbol[hall] = data["hall_symbol"]
                _Symmetry.hall_to_IT_number[hall] = data["number"]

                # Handle Hall to spacegroup using the full H-M name
                key = "international_full"
                name = data[key]

                # Default setting if there are multiple, leave unadorned
                if choice in ("2", "b", "b1", "H"):
                    if (
                        hall in _Symmetry.hall_to_spgname
                        and name != _Symmetry.hall_to_spgname[hall]
                    ):
                        raise RuntimeError(
                            f"{hall=} {key} --> {name} exists: "
                            f"{_Symmetry.hall_to_spgname[hall]}"
                        )
                    else:
                        _Symmetry.hall_to_spgname[hall] = name
                else:
                    # Add other settings to the spacegroup name
                    if choice != "":
                        name += f" :{choice}"
                    if hall not in _Symmetry.hall_to_spgname:
                        # pws if choice != "":
                        # pws print(f"{hall} = {name}   QQQQQQQQQQQQQQQQQQQQQQ")
                        _Symmetry.hall_to_spgname[hall] = name

                # Now handle all the rest
                for key in (
                    "international_full",
                    "international",
                    "international_short",
                    "hall_symbol",
                ):
                    name = data[key]

                    # SPGLib encodes the international name: 'C 2/c = B 2/n 1 1',
                    if key == "international":
                        name = name.split("=")[0].strip()

                    # Default setting if there are multiple, leave unadorned
                    if choice in ("2", "b", "b1", "H"):
                        _Symmetry.spgname_to_hall[name] = hall
                        _Symmetry.spgname_to_system[name] = system_name[key]
                        for txt in ("_", " "):
                            tmp = name.replace(txt, "")
                            _Symmetry.spgname_to_hall[tmp] = hall
                            _Symmetry.spgname_to_system[tmp] = system_name[key]
                            if tmp[-2:] == ":H":
                                tmp = tmp[:-2].strip()
                                _Symmetry.spgname_to_hall[tmp] = hall
                                _Symmetry.spgname_to_system[tmp] = system_name[key]

                    if (
                        key
                        in (
                            "international",
                            "international_short",
                            "international_full",
                        )
                        and choice != ""
                    ):
                        name += f" :{choice}"

                    _Symmetry.spgname_to_hall[name] = hall
                    _Symmetry.spgname_to_system[name] = system_name[key]
                    for txt in ("_", " "):
                        tmp = name.replace(txt, "")
                        _Symmetry.spgname_to_hall[tmp] = hall
                        _Symmetry.spgname_to_system[tmp] = system_name[key]
                        if tmp[-2:] == ":H":
                            tmp = tmp[:-2].strip()
                            _Symmetry.spgname_to_hall[tmp] = hall
                            _Symmetry.spgname_to_system[tmp] = system_name[key]
            # print("spgname_to_hall")
            # pprint.pprint(_Symmetry.spgname_to_hall)
        return _Symmetry.spgname_to_hall

    @staticmethod
    def spacegroup_numbers_to_hall():
        """List of the Hall spacegroup names for the IT number."""
        if _Symmetry.spgno_to_hall is None:
            # Initialize the symmetry data
            _Symmetry.spgno_to_hall = {}
            if _Symmetry.spgname_to_system is None:
                _Symmetry.spgname_to_system = {}
            for hall in range(1, 530):
                data = spglib.get_spacegroup_type(hall)
                spgno = int(data["number"])
                if spgno not in _Symmetry.spgno_to_hall:
                    _Symmetry.spgno_to_hall[spgno] = hall
                    _Symmetry.spgname_to_system[spgno] = "Int_Tables_number"
                    _Symmetry.spgname_to_system[str(spgno)] = "Int_Tables_number"

        return _Symmetry.spgno_to_hall

    @staticmethod
    def to_hall(name):
        """Hall number given full spacegroup name or number."""
        if isinstance(name, int):
            return _Symmetry.spacegroup_numbers_to_hall()[name]
        return _Symmetry.spacegroup_names_to_hall()[name]

    @staticmethod
    def spacegroup_names_to_system():
        """Dictionary of system (name_Hall, name_H-M_full,...) for spacegroup names."""
        if _Symmetry.spgname_to_hall is None:
            _Symmetry.spacegroup_names_to_hall()
        if _Symmetry.spgno_to_hall is None:
            _Symmetry.spacegroup_numbers_to_hall()

        return _Symmetry.spgname_to_system

    @property
    def configuration(self):
        """Return the configuration."""
        return self._configuration

    @property
    def cursor(self):
        return self.system_db.cursor

    @property
    def db(self):
        return self.system_db.db

    @property
    def atom_generators(self):
        """The symmetry operations that create the symmetric atoms."""
        if self._atom_generators is None:
            self._expand()
        return self._atom_generators

    @property
    def atom_to_asymmetric_atom(self):
        """The asymmetric atom related to each atom."""
        if self._atom_to_asymmetric_atom is None:
            self._expand()
        return self._atom_to_asymmetric_atom

    @property
    def bonds_for_asymmetric_bonds(self):
        """List of bonds for each asymmetric bond."""
        if self.n_symops == 1:
            result = [[i] for i in range(self.configuration.bonds.n_asymmetric_bonds)]
        else:
            result = [[] for i in range(self.configuration.bonds.n_asymmetric_bonds)]
            to_asym = self.bond_to_asymmetric_bond
            for i, asym_bond in enumerate(to_asym):
                result[asym_bond].append(i)
        return result

    @property
    def bond_to_asymmetric_bond(self):
        """The asymmetric bond corresponding to each bond."""
        if self._bond_to_asymmetric_bond is None:
            self._expand_bonds()
        return self._bond_to_asymmetric_bond

    @property
    def bond_atoms(self):
        """The pair of symmetric atoms in each bond."""
        if self._bond_atoms is None:
            self._expand_bonds()
        return self._bond_atoms

    @property
    def bond_offsets(self):
        """The triplet of cell offsets of second atom in each bond."""
        if self._bond_offsets is None:
            self._expand_bonds()
        return self._bond_offsets

    @property
    def group(self):
        """The point or space group of the system"""
        self.cursor.execute('SELECT "group" FROM symmetry WHERE id = ?', (self.id,))
        return self.cursor.fetchone()[0]

    @group.setter
    def group(self, value):
        if value == "":
            self.Ws = None
            symops = None
        else:
            if self.configuration.periodicity == 0:
                if value != "C1":
                    raise NotImplementedError("Point groups not implemented yet!")
                symops = ["x,y,z"]
                W4 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

                W4s = []
                W4s.append(W4)
                self._operators = np.array(W4s, dtype=float)
            else:
                # Get the 4x4 augmented matrices
                hall = self.to_hall(value)
                value = self.hall_to_spacegroup_name(hall)
                data = spglib.get_symmetry_from_database(hall)
                W4s = []
                for W, w in zip(
                    data["rotations"].tolist(), data["translations"].tolist()
                ):
                    W4 = []
                    for Wrow, w_i in zip(W, w):
                        W4.append([*Wrow, w_i])
                    W4.append([0, 0, 0, 1])
                    W4s.append(W4)
                self._operators = np.array(W4s, dtype=float)

                # The shorthand strings for the symmetry operations
                symops = []
                for W4 in W4s:
                    symop = []
                    for row in W4[0:3]:
                        line = ""
                        for r, xyz in zip(row[0:3], ("x", "y", "z")):
                            if r == 0:
                                pass
                            elif r == 1:
                                if line != "":
                                    line += "+"
                                line += xyz
                            elif r == -1:
                                line += "-" + xyz
                            else:
                                raise RuntimeError(f"bad rotation: '{row}'")
                        if row[3] == 0:
                            pass
                        else:
                            f = Fraction(row[3]).limit_denominator(10)
                            line += f"{f.numerator:+d}/{f.denominator}"
                        symop.append(line)
                    symops.append(",".join(symop))

        self.symops = symops
        self.db.execute(
            'UPDATE symmetry SET "group" = ? WHERE id = ?', (value, self.id)
        )
        self.db.commit()

    @property
    def hall_number(self):
        if self.configuration.periodicity == 0:
            return ""
        else:
            return self.to_hall(self.group)

    @property
    def hall_symbol(self):
        if self.configuration.periodicity == 0:
            return ""
        else:
            if _Symmetry.hall_to_hall_symbol is None:
                self.spacegroup_names_to_hall
            return _Symmetry.hall_to_hall_symbol[self.to_hall(self.group)]

    @property
    def id(self):
        """The id of this cell."""
        return self._id

    @property
    def inverse_operations(self):
        """The list of inverse symmetry operations for each symop."""
        return [x.index(0) for x in self.symop_products]

    @property
    def IT_number(self):
        if self.configuration.periodicity == 0:
            return ""
        else:
            if _Symmetry.hall_to_hall_symbol is None:
                self.spacegroup_names_to_hall
            return _Symmetry.hall_to_IT_number[self.to_hall(self.group)]

    @property
    def loglevel(self):
        """The logging level for this module."""
        result = logger.getEffectiveLevel()
        tmp = logging.getLevelName(result)
        if "Level" not in tmp:
            result = tmp
        return result

    @loglevel.setter
    def loglevel(self, value):
        logger.setLevel(value)

    @property
    def n_symops(self):
        """The number of symmetry operations."""
        return len(self.symops)

    @property
    def symmetry_matrices(self):
        """The symmetry operations as Numpy matrices (4 x 4 x n_ops)"""
        if self._operators is None:
            # Create the 4x4 matrices
            regexp = re.compile(
                r"([-+]?[0-9]/[0-9])?([-+]?[xyz][-+]?[xyz]?)([-+]?[0-9]/[0-9])?"
            )
            if len(self.symops) == 0:
                raise RuntimeError("The group is not defined!")
            W4s = []
            for line in self.symops:
                W4 = [
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
                for row, op in enumerate(line.split(",")):
                    op = op.strip()
                    match = regexp.fullmatch(op)
                    if match is None:
                        raise ValueError(
                            f"{op} is not a valid symmetry operation ({line})"
                        )
                    m1, m2, m3 = match.groups()
                    for column, xyz in enumerate(("x", "y", "z")):
                        if xyz in m2:
                            if "-" + xyz in m2:
                                W4[row][column] = -1
                            else:
                                W4[row][column] = 1
                    if m1 is None:
                        if m3 is None:
                            pass
                        else:
                            W4[row][3] = float(Fraction(m3))
                    else:
                        W4[row][3] = float(Fraction(m1))
                W4s.append(W4)
            self._operators = np.array(W4s, dtype=float)

        return self._operators

    @property
    def symop_products(self):
        """Matrix of symop ids that are the products of two other symops"""
        if self._symop_products is None:
            self._create_symop_products()
        return self._symop_products

    @property
    def symops(self):
        """The list of shorthand symmetry operations."""
        self.cursor.execute("SELECT symops FROM symmetry WHERE id = ?", (self.id,))
        text = self.cursor.fetchone()[0]

        result = []
        if text != "":
            result = text.split(" | ")
        return result

    @symops.setter
    def symops(self, ops):
        if ops is None:
            text = ""
        else:
            text = " | ".join(ops)
        # Ensure the lower case letters are used! x, y, z not X, Y, Z
        self.db.execute(
            "UPDATE symmetry SET symops = ? WHERE id = ?", (text.lower(), self.id)
        )
        # Unset the group
        self.db.execute('UPDATE symmetry SET "group" = "" WHERE id = ?', (self.id,))
        self.db.commit()

        self._operators = None
        self._symop_products = None
        self._atom_generators = None
        self._atom_to_asymmetric_atom = None
        self._symop_to_atom = None
        self._bond_to_asymmetric_bond = None
        self._bond_atoms = None
        self._bond_offsets = None

    @property
    def symop_to_atom(self):
        """The list of sym ops for each asymmetric atom to create symmetric atoms."""
        if self._symop_to_atom is None:
            self._expand()
        return self._symop_to_atom

    @property
    def system(self):
        """Return the System object that contains this cell."""
        return self._system

    @property
    def system_db(self):
        """Return the SystemDB object that contains this cell."""
        return self._system_db

    def find_spacegroup_from_operators(self):
        """Find the spacegroup from the symmetry operators."""
        if self.configuration.periodicity > 0:
            # Treat P1 as special to avoid bug in spglib...
            if self.n_symops == 1:
                hall_number = 1
                international_number = 1
            else:
                rotations = self.symmetry_matrices[:, 0:3, 0:3]
                translations = self.symmetry_matrices[:, 0:3, 3]
                data = spglib.get_spacegroup_type_from_symmetry(
                    rotations.tolist(),
                    translations.tolist(),
                    self.configuration.cell.vectors(),
                )
                hall_number = data["hall_number"]
                international_number = data["number"]
                while True:
                    data = self.configuration.get_symmetry_data(hall_number)
                    if data is None:
                        raise RuntimeError("Error finding spacegroup from operators.")
                    if data["number"] != international_number:
                        raise RuntimeError(
                            "Error finding setting for spacegroup number "
                            f"{international_number}."
                        )
                    if all(abs(data["origin_shift"]) < 0.001):
                        break
                    hall_number += 1
            return self.hall_to_spacegroup_name(hall_number)

    def reset_atoms(self):
        """The atoms have changed, so need to recalculate the symmetric atoms."""
        self._operators = None
        self._atom_generators = None
        self._atom_to_asymmetric_atom = None
        self._symop_to_atom = None

    def reset_bonds(self):
        """The bonds have changed, so need to recalculate the symmetric bonds."""
        self._bond_to_asymmetric_bond = None
        self._bond_atoms = None
        self._bond_offsets = None

    def symmetrize_atomic_scalar(self, v_sym):
        """Return the symmetrized scalars and deltas.

        Parameters
        ----------
        v_sym : [n_atoms * float] or numpy.ndarray
            The full set of scalar values for the symmetric atoms.

        Returns
        -------
        ([n_asymmetric_atoms * float] or numpy.ndarray, ditto)
            The asymmetric (symmetry unique) values in the same form as given,
            and the delta of the values from the average (or None for C1/P1)
        """
        if self.configuration.periodicity == 0:
            if self.n_symops == 1:
                return v_sym, None
            else:
                raise NotImplementedError("symmetrize_atomic_scalar for molecules!")

        if self.n_symops == 1:
            return v_sym, None

        v_in = np.array(v_sym)
        # print(f"{v_in=}")
        generators = self.atom_generators
        v = np.ndarray((len(generators),), dtype=float)
        delta = np.zeros_like(v_in)
        start = 0
        for asym_atom, ops in enumerate(generators):
            # print(f"{asym_atom=} {ops}")
            n = len(ops)
            # print(f"{n=}")
            v[asym_atom] = np.average(v_in[start : start + n], axis=0)
            delta[start : start + n] = v_in[start : start + n] - v[asym_atom]
            start += n

        if isinstance(v_sym, np.ndarray):
            return v, delta
        else:
            return v.tolist(), delta.tolist()

    def symmetrize_bond_scalar(self, v_sym):
        """Return the symmetrized scalar value and deltas for a bond property

        Parameters
        ----------
        v_sym : [n_bonds * float] or numpy.ndarray
            The full set of scalar values for the symmetric bonds.

        Returns
        -------
        ([n_asymmetric_bonds * float] or numpy.ndarray, ditto)
            The asymmetric (symmetry unique) values in the same form as given,
            and the delta of the values from the average (or None for C1/P1)
        """
        if self.configuration.periodicity == 0:
            if self.n_symops == 1:
                return v_sym, None
            else:
                raise NotImplementedError("symmetrize_bond_scalar for molecules!")

        if self.n_symops == 1:
            return v_sym, None

        v_in = np.array(v_sym)
        generators = self.bonds_for_asymmetric_bonds
        v = np.ndarray((len(generators),), dtype=float)
        delta = np.zeros_like(v_in)
        start = 0
        for asym_bond, ops in enumerate(generators):
            n = len(ops)
            v[asym_bond] = np.average(v_in[start : start + n], axis=0)
            delta[start : start + n] = v_in[start : start + n] - v[asym_bond]
            start += n

        if isinstance(v_sym, np.ndarray):
            return v, delta
        else:
            return v.tolist(), delta.tolist()

    def symmetrize_coordinates(self, xyz_sym, fractionals=True):
        """Return the symmetrized asymmetric coordinates and rms error.

        Parameters
        ----------
        xyz_sym : [n_atoms * [float]] or numpy.ndarray
            The full set of coordinates for the symmetric atoms.
        fractionals : bool = True
            Whether fractional (True) or Cartesian (False) coordinates are given.
            Only important for periodic systems, since molecules use Cartesians.

        Returns
        -------
        ([n_asymmetric_atoms * [float]] or numpy.ndarray, ditto)
            The asymmetric (symmetry unique) coordinates in the same form as given,
            and the delta of the atom positions from the average (or None for C1/P1)
        """
        if self.configuration.periodicity == 0:
            if self.n_symops == 1:
                return xyz_sym, None
            else:
                raise NotImplementedError("symmetrize_coordinates for molecules!")

        if self.n_symops == 1:
            return xyz_sym, None

        # We need fractional coordinates
        if fractionals:
            uvw_sym = np.array(xyz_sym)
        else:
            uvw_sym = self.configuration.cell.to_fractionals(xyz_sym, as_array=True)

        # Translate into cell
        uvw_sym -= np.floor(uvw_sym)

        generators = self.atom_generators
        inverse_ops = self.inverse_operations
        uvw = np.ndarray((len(generators), 3), dtype=float)
        delta = np.zeros_like(uvw_sym)
        atom = 0
        for asym_atom, ops in enumerate(generators):
            tmp = np.ndarray((len(ops), 3), dtype=float)
            start = atom
            for i, op in enumerate(ops):
                tmp[i] = self.vector_x_symop(uvw_sym[atom], inverse_ops[op])
                atom += 1
            tmp -= np.floor(tmp)
            uvw[asym_atom] = np.average(tmp, axis=0)
            delta[start:atom] = tmp - uvw[asym_atom]

        if isinstance(xyz_sym, np.ndarray):
            if fractionals:
                return uvw.round(8), delta.round(8)
            else:
                tmp = self.configuration.cell.to_cartesians(uvw, as_array=True)
                return tmp.round(8), delta.round(8)
        else:
            if fractionals:
                return uvw.round(8).tolist(), delta.round(8).tolist()
            else:
                tmp = self.configuration.cell.to_cartesians(uvw, as_array=True)
                return tmp.round(8).tolist(), delta.round(8).tolist()

    def update_group(self, value):
        if self.configuration.periodicity == 0:
            if value != "C1":
                raise NotImplementedError("Point groups not implemented yet!")
        self.db.execute(
            'UPDATE symmetry SET "group" = ? WHERE id = ?', (value, self.id)
        )
        self.db.commit()

    def vector_x_symop(self, vector, symop, translation=True):
        """Multiply a vector by a symmetry matrix."""
        sym_mat = self.symmetry_matrices[symop]
        if translation:
            v = np.array([0.0, 0.0, 0.0, 1.0])
        else:
            v = np.array([0.0, 0.0, 0.0, 0.0])
        v[0:3] = vector
        xformed = np.einsum("ij,j", sym_mat, v)
        return xformed[0:3]

    def _expand(self):
        """Setup the information for going from asymmetric cell to full cell.

        There is an issue with some CIF files, notably those from the Cambridge
        structural database, where there are atoms in the asymmetric unit that
        are duplicates. It appears that the CIF file is written with entire molecules
        in some cases. So if the molecule is on a symmetry element, the atoms that are
        related by symmetry are explicitly in the file.

        To handle this we need to check for duplicates and remove them. Not a bad check
        to have, anyway.
        """
        operators = self.symmetry_matrices

        symbols = self.configuration.atoms.asymmetric_symbols
        n_atoms = self.configuration.n_asymmetric_atoms
        logger.info(f"Expanding {n_atoms} asymmetric atoms to full cell.")
        if n_atoms == 0:
            self._atom_generators = None
            self._symop_to_atom = None
        self._atom_generators = []  # Symmetry ops that create the symmetric atoms
        self._symop_to_atom = []  # The symmetry atom resulting from symmetry ops
        if self.configuration.periodicity == 3:
            if len(operators) == 1:
                # P1 is a special case. Nothing to do.
                self._atom_generators = [[0] for i in range(n_atoms)]
                self._symop_to_atom = [[i] for i in range(n_atoms)]
                self._atom_to_asymmetric_atom = [*range(n_atoms)]
            else:
                uvw0 = self.configuration.atoms.get_coordinates(
                    as_array=True, asymmetric=True
                )
                if uvw0.shape[0] != n_atoms:
                    raise RuntimeError(
                        f"Mismatch of number of atoms in symmetry: {uvw0.shape[0]} != "
                        f"{n_atoms}"
                    )
                logger.debug(
                    f"""
{n_atoms=}
{uvw0.shape=}"
Original coordinates
{uvw0}
                    """
                )

                uvw = np.ndarray((n_atoms, 4))
                uvw[:, 0:3] = uvw0[:, :]
                logger.debug(f"\nExpanded coordinates\n{uvw}")
                uvw[:, 3] = 1
                # logger.debug(self.symops)
                logger.debug(f"\n{operators.shape=}")
                # logger.debug(operators)
                # logger.debug(f"{uvw.shape=}")
                # logger.debug("Expanded coordinates")
                # logger.debug(uvw)
                xformed = np.einsum("ijk,lk", operators, uvw)
                logger.debug(f"\n{xformed.shape=}\n{xformed}")

                # For comparison, bring all atoms into the cell [0..1)
                tmp = xformed - np.floor(xformed)

                logger.debug(f"\ntmp = Transformed and floored coordinates\n{tmp}")

                # Keep track of atoms and coordinates found so can check for duplicates
                found = {}

                n_sym_atoms = 0
                n_asymmetric_atoms = 0
                for i in range(n_atoms):
                    values, I1, I2 = np.unique(
                        np.round(tmp[:, :3, i], 4),
                        axis=0,
                        return_index=True,
                        return_inverse=True,
                    )
                    # pprint.pprint(np.round(tmp[:, :3, i], 4).tolist())
                    logger.debug(
                        f"""
{i=}: {symbols[i]}
{np.round(tmp[:, :3, i], 4)}
nvalues
{values}
I1
{I1}
I2
{I2}
                        """
                    )
                    # Check for duplicates
                    if tuple(values[0]) in found:
                        if found[tuple(values[0])] != symbols[i]:
                            raise RuntimeError(
                                "Duplicate atoms with different symbols: "
                                f"{symbols[i]} != {found[tuple(values[0])]} at "
                                f"{values[0]}"
                            )
                        else:
                            # Same element, so ignore
                            logger.debug("\nDuplicate atom, ignoring")
                            continue

                    for value in values:
                        found[tuple(value)] = symbols[i]

                    n_asymmetric_atoms += 1
                    I2 += n_sym_atoms
                    self._atom_generators.append(I1.tolist())
                    self._symop_to_atom.append(I2.tolist())
                    n_sym_atoms += I1.shape[0]
                tmp = []
                for i, generators in enumerate(self._atom_generators):
                    tmp.extend([i] * len(generators))
                self._atom_to_asymmetric_atom = tmp
                logger.info(f"{n_asymmetric_atoms=} {n_sym_atoms=}")

    def _expand_bonds(self):
        """Expand the list of asymmetric bonds to the full list.

        Keep track of the following:

            1. bond_to_asymmetric_bond: the asymmetric bond index for each bond
            2. bond_atoms: The pairs of symmetric atoms that form the bond
            3. bond_offsets: Triplets of offsets for bonds
        """
        logger.debug("In expand_bonds")

        if self._bond_to_asymmetric_bond is not None:
            return

        self._bond_to_asymmetric_bond = []
        self._bond_atoms = []
        self._bond_offsets = []

        bonds = self.configuration.bonds
        if bonds.n_asymmetric_bonds == 0:
            return

        indx = {j: i for i, j in enumerate(self.configuration.atoms.ids)}

        Is = [indx[i] for i in bonds.get_column_data("i")]
        Js = [indx[j] for j in bonds.get_column_data("j")]
        symop1s = bonds.get_column_data("symop1")
        symop2s = bonds.get_column_data("symop2")

        symop_to_atom = self.symop_to_atom
        product = self.symop_products
        n_symops = self.n_symops

        # Need to coordinates to sort out offsets introduced by symmetry
        uvw = self.configuration.atoms.get_coordinates(as_array=True)
        uvw_sym = self.configuration.atoms.get_coordinates(
            as_array=True, asymmetric=True
        )

        logger.debug(f"Asymmetric Coordinates:\n{str(uvw_sym)}\n\n")
        logger.debug(f"Coordinates:\n{str(uvw)}\n\n")

        asym_bonds = {}
        found = []
        # First work out offsets and look for bonds sharing atoms but different offsets
        for i, j, symop1, symop2 in zip(Is, Js, symop1s, symop2s):
            if symop1 == ".":
                op1 = 0
                offset1 = [0, 0, 0]
            else:
                if "_" in symop1:
                    op1, tmp = symop1.split("_")
                    op1 = int(op1) - 1
                    offset1 = [int(q) - 5 for q in tmp]
                else:
                    op1 = int(symop1) - 1
                    offset1 = [0, 0, 0]
            if symop2 == ".":
                op2 = 0
                offset2 = [0, 0, 0]
            else:
                if "_" in symop2:
                    op2, tmp = symop2.split("_")
                    op2 = int(op2) - 1
                    offset2 = [int(q) - 5 for q in tmp]
                else:
                    op2 = int(symop2) - 1
                    offset2 = [0, 0, 0]

            # Check if we've seen this pair, but with different offset.
            if i < j:
                key = (i, op1, j, op2)
                offkey = (offset1, offset2)
            else:
                key = (j, op2, i, op1)
                offkey = (offset2, offset1)
            logger.debug(f"{key=} {offkey}")
            if key in asym_bonds:
                if offkey in asym_bonds[key]:
                    raise RuntimeError(
                        f"Duplicate bonds specified for {i} ({symop1}) --  "
                        f"{j} ({symop2})"
                    )
                else:
                    asym_bonds[key].append(offkey)
            else:
                asym_bonds[key] = [offkey]

        bond_no = -1
        for ij, offsets in asym_bonds.items():
            i, op1, j, op2 = ij
            use_offsets = len(offsets) > 1 or (i == j and op1 == op2)
            for tmp in offsets:
                bond_no += 1
                offset1, offset2 = tmp
                if use_offsets:
                    logger.debug(
                        f"bond {bond_no}: {i}({op1}) - {j}({op2}) {offset1} {offset2}"
                    )
                else:
                    logger.debug(
                        f"bond {bond_no}: {i}({op1}) - {j}({op2}) ignore offsets"
                    )

                for op in range(n_symops):
                    prod_op1 = product[op1][op]
                    prod_op2 = product[op2][op]
                    iatom = symop_to_atom[i][prod_op1]
                    jatom = symop_to_atom[j][prod_op2]
                    ioff = self.vector_x_symop(offset1, prod_op1, translation=False)
                    joff = self.vector_x_symop(offset2, prod_op2, translation=False)

                    if iatom == jatom:
                        logger.debug(
                            f"iatom == jatom ({iatom}, {jatom}) {use_offsets=}"
                        )
                        if use_offsets:
                            if np.all(offset1 == offset2):
                                logger.debug(f"{offset1=} {offset2=}")
                                continue
                        else:
                            continue

                    if iatom > jatom:
                        ii, jj = j, i
                        iatom, jatom = jatom, iatom
                        offset1, offset2 = offset2, offset1
                        prod_op1, prod_op2 = prod_op2, prod_op1
                    else:
                        ii, jj = i, j

                    uvw1 = self.vector_x_symop(uvw_sym[ii], prod_op1)
                    uvw2 = self.vector_x_symop(uvw_sym[jj], prod_op2)

                    delta = uvw2 - uvw1
                    delta = delta.round(4)
                    # Not sure is <= or < :-)
                    tmpoff = np.select([delta < -0.5, delta > 0.5], [1, -1], 0)
                    if use_offsets:
                        ioff = np.array(offset1)
                        joff = np.array(tmpoff) + offset2
                    else:
                        ioff = np.array((0, 0, 0))
                        joff = np.array(tmpoff)

                    delta += joff - ioff

                    if np.any(abs(delta) > 1.0):
                        logger.debug(f"Throwing out delta = {str(delta)}")
                        continue

                    if (iatom, jatom, delta.tolist()) in found:
                        continue

                    found.append((iatom, jatom, delta.tolist()))

                    self._bond_to_asymmetric_bond.append(bond_no)
                    self._bond_atoms.append((iatom, jatom))
                    self._bond_offsets.append((ioff.tolist(), joff.tolist()))

                    if logger.isEnabledFor(logging.DEBUG):
                        uvw1 = self.vector_x_symop(uvw_sym[ii], prod_op1)
                        uvw2 = self.vector_x_symop(uvw_sym[jj], prod_op2)

                        delta = uvw2 - uvw1
                        tmpoff = np.select([delta <= -0.5, delta > 0.5], [1, -1], 0)
                        delta += joff - ioff
                        delta_xyz = self.configuration.cell.to_cartesians(
                            delta, as_array=True
                        )
                        r = np.round(np.linalg.norm(delta_xyz), 4)
                        xi, yi, zi = uvw1.tolist()
                        xj, yj, zj = uvw2.tolist()
                        xd, yd, zd = delta.tolist()
                        off1x, off1y, off1z = ioff.tolist()
                        off2x, off2y, off2z = joff.tolist()
                        logger.debug(
                            f"\t{iatom:3} {jatom:3} | {ii} {product[op1][op]:3} "
                            f"{jj} {product[op2][op]:3} | "
                            f"{xi:5.2f} {yi:5.2f} {zi:5.2f}"
                            f"({off1x:4} {off1y:4} {off1z:4}) "
                            f"{xj:5.2f} {yj:5.2f} {zj:5.2f}"
                            f"({off2x:4} {off2y:4} {off2z:4}) "
                            f"| {xd:5.2f} {yd:5.2f} {zd:5.2f}"
                            f" |  {r=:7.4f}"
                        )

    def _create_symop_products(self):
        """Create the matrix of the product of symmetry operators."""
        products = []
        ops = self.symmetry_matrices
        n = self.n_symops
        # Take the i'th symop and multiply by all the rest
        for i in range(n):
            prod = np.einsum("jk,ikl", ops[i, :, :], ops)
            row = []
            # For the i-j product, find it in the original list
            for j in range(n):
                # logger.debug(f"{j}:")

                # Subtract the product from the original symops
                tmp = ops - prod[j, :, :]

                # Shift any translation into the cell, i.e. [0,1)
                tmp[:, :, 3] = tmp[:, :, 3] - np.floor(tmp[:, :, 3])

                # We are looking for zeros, so sum the abs values
                tmp = tmp.reshape(n, 16)
                tmp2 = np.sum(np.abs(tmp), axis=1)

                # And find any matrices that are all zero
                hits = np.where(tmp2 == 0.0)[0]

                if len(hits) != 1:
                    symops = self.symops
                    logger.warning("Symmetry matrices")
                    for k in range(n):
                        logger.warning(f"{k}: {symops[k]}\n{str(ops[k])}")
                    logger.warning("")
                    logger.warning(f"Symop {i}:\n{str(ops[i])}")
                    logger.warning("")
                    logger.warning("Products")
                    for k in range(n):
                        logger.warning(f"\n{str(tmp[k])}")
                    raise RuntimeError(f"Problem with products of symops: {len(hits)=}")
                k = hits[0]
                row.append(k)
                # logger.debug(prod[j])
                # logger.debug(f"symop {k} is")
                # logger.debug(ops[k])
            if len(row) != n:
                raise RuntimeError(f"Not enough products of the symop {i}")
            # logger.debug(f"{i}: ({len(row)}) {row}")
            products.append(row)
        self._symop_products = products

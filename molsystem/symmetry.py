# -*- coding: utf-8 -*-

import json  # noqa: F401
import logging

import numpy as np
import spglib

logger = logging.getLogger(__name__)
# logger.setLevel("DEBUG")


class _Symmetry(object):
    """A class to handle point and space group symmetry."""

    spgno_to_hall = None
    spgname_to_hall = None
    spgname_to_system = None

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
        self._atom_generators = None  # Sym operations that create the symmetric atoms
        self._symop_to_atom = None  # The symmetry atom resulting from symmetry ops
        self._atom_to_asymmetric_atom = None  # The asymetric atom for each atom
        super().__init__()

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
                        elif row[3] == 0.5:
                            line += "+1/2"
                        elif row[3] == -0.5:
                            line += "-1/2"
                        elif row[3] == 0.25:
                            line += "+1/4"
                        elif row[3] == -0.25:
                            line += "-1/4"
                        else:
                            raise RuntimeError(f"bad translation: '{row[3]}'")
                        symop.append(line)
                    symops.append(",".join(symop))

        self.symops = symops
        self.db.execute(
            'UPDATE symmetry SET "group" = ? WHERE id = ?', (value, self.id)
        )
        self.db.commit()

    @property
    def id(self):
        """The id of this cell."""
        return self._id

    @property
    def n_symops(self):
        """The number of symmetry operations."""
        return len(self.symops)

    @property
    def symmetry_matrices(self):
        """The symmetry operations as Numpy matrices (4 x 4 x n_ops)"""
        if self._operators is None:
            # Create the 4x4 matrices
            if len(self.symops) == 0:
                raise RuntimeError("The group is not defined!")
            W4s = []
            for line in self.symops:
                W4 = []
                for op in line.split(","):
                    row = [0, 0, 0, 0]
                    for index, xyz in enumerate(("x", "y", "z")):
                        if "-" + xyz in op:
                            row[index] = -1
                        elif xyz in op:
                            row[index] = 1
                    for text, value in zip(
                        ("-1/2", "1/2", "-1/4", "1/4"), (-0.5, 0.5, -0.25, 0.25)
                    ):
                        if text in op:
                            row[3] = value
                            break
                    W4.append(row)
                W4.append([0, 0, 0, 1])
                W4s.append(W4)
            self._operators = np.array(W4s, dtype=float)
        return self._operators

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
            self._operators = None
        else:
            text = " | ".join(ops)

        self.db.execute("UPDATE symmetry SET symops = ? WHERE id = ?", (text, self.id))
        # Unset the group
        self.db.execute('UPDATE symmetry SET "group" = "" WHERE id = ?', (self.id,))
        self.db.commit()

        self._operators = None
        self._atom_generators = None
        self._atom_to_asymmetric_atom = None
        self._symop_to_atom = None

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

    def to_hall(self, name):
        """Hall number given full spacegroup name or number."""
        if isinstance(name, int):
            return self.spacegroup_numbers_to_hall[name]
        return self.spacegroup_names_to_hall[name]

    @property
    def spacegroup_numbers_to_hall(self):
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

    @property
    def spacegroup_names_to_hall(self):
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
            if _Symmetry.spgname_to_system is None:
                _Symmetry.spgname_to_system = {}
            for hall in range(1, 530):
                data = spglib.get_spacegroup_type(hall)
                choice = data["choice"]
                for key in (
                    "international_full",
                    "international",
                    "international_short",
                    "hall_symbol",
                ):
                    name = data[key]
                    if "international" in key and choice in ("2",):
                        if (
                            name in _Symmetry.spgname_to_hall
                            and hall != _Symmetry.spgname_to_hall[name]
                        ):
                            raise RuntimeError(
                                f"{hall=} {key} --> {name} exists: "
                                f"{_Symmetry.spgname_to_hall[name]}"
                            )
                        _Symmetry.spgname_to_hall[name] = hall
                        _Symmetry.spgname_to_system[name] = system_name[key]
                        name = name.replace("_", "")
                        _Symmetry.spgname_to_hall[name] = hall
                    if choice != "":
                        name += f":{choice}"
                    if (
                        name in _Symmetry.spgname_to_hall
                        and hall != _Symmetry.spgname_to_hall[name]
                    ):
                        raise RuntimeError(
                            f"{hall=} {key} --> {name} exists: "
                            f"{_Symmetry.spgname_to_hall[name]}"
                        )
                    _Symmetry.spgname_to_hall[name] = hall
                    _Symmetry.spgname_to_system[name] = system_name[key]
                    name = name.replace("_", "")
                    _Symmetry.spgname_to_hall[name] = hall
        return _Symmetry.spgname_to_hall

    @property
    def spacegroup_names_to_system(self):
        """Dictionary of system (name_Hall, name_H-M_full,...) for spacegroup names."""
        if _Symmetry.spgname_to_hall is None:
            self.spacegroup_names_to_hall
        if _Symmetry.spgno_to_hall is None:
            self.spacegroup_numbers_to_hall

        return _Symmetry.spgname_to_system

    def _expand(self):
        """Setup the information for going from asymmetric cell to full cell."""
        operators = self.symmetry_matrices

        n_atoms = self.configuration.n_asymmetric_atoms
        if n_atoms == 0:
            self._atom_generators = None
            self._symop_to_atom = None
        self._atom_generators = []  # Symmetry ops that create the symmetric atoms
        self._symop_to_atom = []  # The symmetry atom resulting from symmetry ops
        if self.configuration.periodicity == 3:
            if len(operators) == 1:
                # P1 is a special case. Nothing to do!
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
                logger.debug("Original coordinates")
                logger.debug(uvw0)

                uvw = np.ndarray((n_atoms, 4))
                uvw[:, 0:3] = uvw0[:, :]
                logger.debug("Expanded coordinates")
                logger.debug(uvw)
                uvw[:, 3] = 1
                # logger.debug(self.symops)
                logger.debug(f"{operators.shape=}")
                # logger.debug(operators)
                logger.debug(f"{uvw.shape=}")
                logger.debug("Expanded coordinates")
                logger.debug(uvw)
                xformed = np.einsum("ijk,lk", operators, uvw)
                logger.debug(f"{xformed.shape=}")
                # logger.debug(xformed)

                # For comparison, bring all atoms into the cell [0..1)
                tmp = xformed - np.floor(xformed)

                logger.debug("tmp")
                logger.debug(tmp)

                for i in range(n_atoms):
                    values, I1, I2 = np.unique(
                        np.round(tmp[:, :3, i], 4),
                        axis=0,
                        return_index=True,
                        return_inverse=True,
                    )
                    logger.debug(i)
                    logger.debug("")
                    logger.debug(values)
                    logger.debug("")
                    logger.debug(I1)
                    logger.debug("")
                    logger.debug(I2)
                    self._atom_generators.append(I1)
                    self._symop_to_atom.append(I2)
                tmp = []
                for i, generators in enumerate(self._atom_generators):
                    tmp.extend([i] * len(generators))
                self._atom_to_asymmetric_atom = tmp
        else:
            if len(operators) == 1:
                # P1 is a special case. Nothing to do!
                self._atom_generators = [[0] for i in range(n_atoms)]
                self._symop_to_atom = [[i] for i in range(n_atoms)]
                self._atom_to_asymmetric_atom = [*range(n_atoms)]

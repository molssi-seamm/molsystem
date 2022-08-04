# -*- coding: utf-8 -*-

import logging

import spglib

logger = logging.getLogger(__name__)


class _Symmetry(object):
    """A class to handle point and space group symmetry.

    :meta public:
    """

    spgno_to_hall = None
    spgname_to_hall = None

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
    def group(self):
        """The point or space group of the system"""
        self.cursor.execute('SELECT "group" FROM symmetry WHERE id = ?', (self.id,))
        return self.cursor.fetchone()[0]

    @group.setter
    def group(self, value):
        self.db.execute(
            'UPDATE symmetry SET "group" = ? WHERE id = ?', (value, self.id)
        )
        self.db.commit()

    @property
    def id(self):
        """The id of this cell."""
        return self._id

    @property
    def system(self):
        """Return the System object that contains this cell."""
        return self._system

    @property
    def system_db(self):
        """Return the SystemDB object that contains this cell."""
        return self._system_db

    @classmethod
    def full_spgname_to_hall(cls, spgname, as_strings=False):
        """Hall number given full spacegroup name."""
        if _Symmetry.spgno_to_hall is None:
            # Initialize the symmetry data
            _Symmetry.spgno_to_hall = {}
            _Symmetry.spgname_to_hall = {}
            for hall in range(1, 530):
                data = spglib.get_spacegroup_type(hall)
                spgno = data["number"]
                if spgno not in _Symmetry.spgno_to_hall:
                    _Symmetry.spgno_to_hall[spgno] = hall
                name = data["international_full"]
                if name not in _Symmetry.spgname_to_hall:
                    _Symmetry.spgname_to_hall[name] = hall
                    name = name.replace("_", "")
                    _Symmetry.spgname_to_hall[name] = hall

        return _Symmetry.spgname_to_hall[spgname]

    @classmethod
    def symops_as_strings(cls, spgname):
        """Return the symmetry operators for the group.

        Parameters
        ----------
        spgname : str
            The full International spacegroup name

        Result
        ------
        str
            The spacegroup operators as strings, e.g. 'x, y+1/2, z'
        """
        hall = _Symmetry.full_spgname_to_hall(spgname)
        data = spglib.get_symmetry_from_database(hall)
        result = []
        for rotation, translation in zip(
            data["rotations"].tolist(), data["translations"].tolist()
        ):
            symops = []
            for r1, t in zip(rotation, translation):
                line = ""
                for r, xyz in zip(r1, ("x", "y", "z")):
                    if r == 0:
                        pass
                    elif r == 1:
                        line += xyz
                    elif r == -1:
                        line += "-" + xyz
                    else:
                        raise RuntimeError(f"bad rotation: '{r1}'")

                if t == 0:
                    pass
                elif t == 0.5:
                    line += "+1/2"
                elif t == -0.5:
                    line += "-1/2"
                elif t == 0.25:
                    line += "+1/4"
                elif t == -0.25:
                    line += "-1/4"
                else:
                    raise RuntimeError(f"bad translation: '{translation}'")
                symops.append(line)
            result.append(",".join(symops))

        return result

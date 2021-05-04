# -*- coding: utf-8 -*-

import logging
import math

import numpy

logger = logging.getLogger(__name__)


def cos(value):
    return math.cos(math.radians(value))


def sin(value):
    return math.sin(math.radians(value))


class Cell(object):
    """A class to handle cell parameters and their transformations."""

    def __init__(self, a, b, c, alpha, beta, gamma):
        self._parameters = [a, b, c, alpha, beta, gamma]

    def __getitem__(self, key):
        """Allow [] to access the data!"""
        return self._parameters[key]

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        self._parameters[key] = value

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self._parameters)

    def __len__(self) -> int:
        """The len() command"""
        return len(self._parameters)

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        return self.equal(other, tol=1.0e-12)

    def __repr__(self):
        """The representation of this object"""
        return repr(self._parameters)

    def __str__(self):
        """The pretty string representation of this object"""
        return str(self._parameters)

    @property
    def a(self):
        """The length of the first cell vector."""
        return self._parameters[0]

    @a.setter
    def a(self, value):
        self._parameters[0] = value
        return list(self._parameters)

    @property
    def b(self):
        """The length of the second cell vector."""
        return self._parameters[1]

    @b.setter
    def b(self, value):
        self._parameters[1] = value
        return list(self._parameters)

    @property
    def c(self):
        """The length of the third cell vector."""
        return self._parameters[2]

    @c.setter
    def c(self, value):
        self._parameters[2] = value
        return list(self._parameters)

    @property
    def alpha(self):
        """The angle between b and c."""
        return self._parameters[3]

    @alpha.setter
    def alpha(self, value):
        self._parameters[3] = value
        return list(self._parameters)

    @property
    def beta(self):
        """The angle between a and c."""
        return self._parameters[4]

    @beta.setter
    def beta(self, value):
        self._parameters[4] = value
        return list(self._parameters)

    @property
    def gamma(self):
        """The angle between a and b."""
        return self._parameters[5]

    @gamma.setter
    def gamma(self, value):
        self._parameters[5] = value
        return list(self._parameters)

    @property
    def parameters(self):
        """The cell parameters as a list."""
        return list(self._parameters)

    @parameters.setter
    def parameters(self, value):
        if len(value) != 6:
            raise ValueError("parameters must be of length 6")
        self._parameters = list(value)

    @property
    def volume(self):
        """The volume of the cell."""
        a, b, c, alpha, beta, gamma = self.parameters
        # Roundoff errors!
        value = cos(alpha) * cos(beta) * cos(gamma)
        if value < 0.0 and abs(value) < 1.0e-8:
            value = 0.0
        return a * b * c * (
            1 - cos(alpha) ** 2 - cos(beta) ** 2 - cos(gamma) ** 2
        ) + 2 * math.sqrt(value)

    def equal(self, other, tol=1.0e-6):
        """Check if we are equal to another iterable to within a tolerance.

        Parameters
        ----------
        other : iterable
            The other object to check against
        tol : float = 1.0e-06
            The tolerance for comparing floating point numbers.

        Returns
        -------
        equals : bool
            Boolean indicating whether the two are equal.
        """
        if len(other) != 6:
            return False

        for i, j in zip(self, other):
            if abs(i - j) > tol:
                return False

        return True

    def to_cartesians(self, uvw, as_array=False):
        """Convert fraction coordinates to Cartesians

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        uvw : [N][3*float] or ndarray
            The fractional coordinates.

        Returns
        -------
        xyz : [N][float*3] or ndarray
            The Cartesian coordinates.
        """
        if isinstance(uvw, numpy.ndarray):
            UVW = uvw
        else:
            UVW = numpy.array(uvw)

        T = self.to_cartesians_transform(as_array=True)
        XYZ = UVW @ T

        if as_array:
            return XYZ
        else:
            return XYZ.tolist()

    def to_cartesians_transform(self, as_array=False):
        """Matrix to convert fractional coordinates to Cartesian.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        as_array : bool = False
            Whether to return a numpy array or Python lists

        Returns
        -------
        transform : [N][float*3] or ndarray
            The transformation matrix
        """
        a, b, c, alpha, beta, gamma = self.parameters

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sg = sin(gamma)

        V = a * b * c * math.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg)
        # Transpose of ...
        # [a, b * cg, c * cb],
        # [0, b * sg, c * (ca - cb * cg) / sg],
        # [0, 0, V / (a * b * sg)]
        T = [
            [a, 0, 0],
            [b * cg, b * sg, 0],
            [c * cb, c * (ca - cb * cg) / sg, V / (a * b * sg)],
        ]  # yapf: disable

        if as_array:
            return numpy.array(T)
        else:
            return T

    def to_fractionals(self, xyz, as_array=False):
        """Convert Cartesian coordinates to fractional.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        xyz : [N][3*float] or ndarray
            The Cartesian coordinates.

        Returns
        -------
        uvw : [N][float*3] or ndarray
            The ractional coordinates.
        """
        if isinstance(xyz, numpy.ndarray):
            XYZ = xyz
        else:
            XYZ = numpy.array(xyz)

        T = self.to_fractionals_transform(as_array=True)
        UVW = XYZ @ T

        if as_array:
            return UVW
        else:
            return UVW.tolist()

    def to_fractionals_transform(self, as_array=False):
        """Matrix to convert Cartesian coordinates to fractional.

        see https://en.wikipedia.org/wiki/Fractional_coordinates for a
        description.

        Parameters
        ----------
        as_array : bool = False
            Whether to return a numpy array or Python lists

        Returns
        -------
        transform : [N][float*3] or ndarray
            The transformation matrix
        """
        a, b, c, alpha, beta, gamma = self.parameters

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sg = sin(gamma)

        V = a * b * c * math.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg)
        # Transpose...
        # [1 / a, -cg / (a * sg), b * c * (ca * cg - cb) / (V * sg)],
        # [0, 1 / (b * sg), a * c * (cb * cg - ca) / (V * sg)],
        # [0, 0, a * b * sg / V]
        T = [
            [1 / a, 0, 0],
            [-cg / (a * sg), 1 / (b * sg), 0],
            [
                b * c * (ca * cg - cb) / (V * sg),
                a * c * (cb * cg - ca) / (V * sg),
                a * b * sg / V,
            ],
        ]  # yapf: disable

        if as_array:
            return numpy.array(T)
        else:
            return T


class _Cell(Cell):
    """A class for handling cell parameters as part of MolSystem.

    Provides all the functionality of the Cell class, but keeps the cell
    data in the database.

    :meta public:
    """

    def __init__(self, configuration):
        """Initialize from the database.

        Parameters
        ----------
        system_db : SystemDB
            The SystemDB instance that we are working with.
        _id : int
            The id of this particular cell.
        """
        self._configuration = configuration
        self._system = self._configuration.system
        self._system_db = self._system.system_db
        self._id = configuration.cell_id

        self.cursor.execute(
            "SELECT a, b, c, alpha, beta, gamma FROM cell WHERE id = ?", (self._id,)
        )
        super().__init__(*self.cursor.fetchone())

    def __enter__(self):
        """Copy the tables to a backup for a 'with' statement."""
        self.system_db["cell"].__enter__()
        return self

    def __exit__(self, etype, value, traceback):
        """Handle returning from a 'with' statement."""
        if etype is None:
            self.configuration.version = self.configuration.version + 1
        return self.system_db["cell"].__exit__(etype, value, traceback)

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        # This gets rid if LGTM warning...
        return self.equal(other, tol=1.0e-12)

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        self._parameters[key] = value
        self._save()

    @property
    def cursor(self):
        return self.system_db.cursor

    @property
    def db(self):
        return self.system_db.db

    @property
    def configuration(self):
        """Return the configuration."""
        return self._configuration

    @property
    def id(self):
        """The id of this cell."""
        return self._id

    @property
    def a(self):
        """The length of the first cell vector."""
        return self._parameters[0]

    @a.setter
    def a(self, value):
        self._parameters[0] = value
        self._save()
        return list(self._parameters)

    @property
    def b(self):
        """The length of the second cell vector."""
        return self._parameters[1]

    @b.setter
    def b(self, value):
        self._parameters[1] = value
        self._save()
        return list(self._parameters)

    @property
    def c(self):
        """The length of the third cell vector."""
        return self._parameters[2]

    @c.setter
    def c(self, value):
        self._parameters[2] = value
        self._save()
        return list(self._parameters)

    @property
    def alpha(self):
        """The angle between b and c."""
        return self._parameters[3]

    @alpha.setter
    def alpha(self, value):
        self._parameters[3] = value
        self._save()
        return list(self._parameters)

    @property
    def beta(self):
        """The angle between a and c."""
        return self._parameters[4]

    @beta.setter
    def beta(self, value):
        self._parameters[4] = value
        self._save()
        return list(self._parameters)

    @property
    def gamma(self):
        """The angle between a and b."""
        return self._parameters[5]

    @gamma.setter
    def gamma(self, value):
        self._parameters[5] = value
        self._save()
        return list(self._parameters)

    @property
    def parameters(self):
        """The cell parameters as a list."""
        return list(self._parameters)

    @parameters.setter
    def parameters(self, value):
        if len(value) != 6:
            raise ValueError("parameters must be of length 6")
        self._parameters = list(value)
        self._save()
        return list(self._parameters)

    @property
    def system(self):
        """Return the System object that contains this cell."""
        return self._system

    @property
    def system_db(self):
        """Return the SystemDB object that contains this cell."""
        return self._system_db

    def diff(self, other):
        """Difference between these cell and another

        Parameters
        ----------
        other : _Cell
            The other cell to diff against

        Result
        ------
        result : Dict
            The differences, described in a dictionary
        """
        result = {}

        # Check the columns
        columns = self._columns()
        other_columns = other._columns()

        column_defs = ", ".join(columns)
        other_column_defs = ", ".join(other_columns)

        if columns == other_columns:
            column_def = column_defs
        else:
            added = columns - other_columns
            if len(added) > 0:
                result["columns added"] = list(added)
            deleted = other_columns - columns
            if len(deleted) > 0:
                result["columns deleted"] = list(deleted)

            in_common = other_columns & columns
            if len(in_common) > 0:
                column_def = ", ".join(in_common)
            else:
                # No columns shared
                return result

        # Need to check the contents of the tables. See if they are in the same
        # database or if we need to attach the other database temporarily.
        db = self.system_db
        other_db = other.system_db

        detach = False
        schema = self.schema
        if db.filename != other_db.filename:
            if db.is_attached(other_db):
                other_schema = db.attached_as(other_db)
            else:
                # Attach the other system_db in order to do comparisons.
                other_schema = self.system_db.attach(other_db)
                detach = True
        else:
            other_schema = other.schema

        _id = self.id
        other_id = other.id

        changed = {}
        last = None
        sql = f"""
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {other_schema}.bond
           WHERE id = {other_id}
          EXCEPT
          SELECT {column_def}
            FROM {schema}.bond
           WHERE id = {_id}
        )
         UNION ALL
        SELECT * FROM
        (
          SELECT {column_def}
            FROM {schema}.bond
           WHERE id = {_id}
          EXCEPT
          SELECT {column_def}
            FROM {other_schema}.bond
           WHERE id = {other_id}
        )
        ORDER BY id
        """

        for row in self.db.execute(sql):
            if last is None:
                last = row
            elif row["id"] == last["id"]:
                # changes = []
                changes = set()
                for k1, v1, v2 in zip(last.keys(), last, row):
                    if v1 != v2:
                        changes.add((k1, v1, v2))
                changed[row["id"]] = changes
                last = None
            else:
                last = row
        if len(changed) > 0:
            result["changed"] = changed

        # See about the rows added
        added = {}
        sql = f"""
        SELECT {column_defs}
          FROM {schema}.bond
         WHERE id = {_id}
           AND id <> {other_id}
        """
        for row in self.db.execute(sql):
            added[row["id"]] = row[1:]

        if len(added) > 0:
            result["columns in added rows"] = row.keys()[1:]
            result["added"] = added

        # See about the rows deleted
        deleted = {}
        sql = f"""
        SELECT {other_column_defs}
          FROM {other_schema}.bond
         WHERE id = {other_id}
           AND id <> {_id}
        """
        for row in self.db.execute(sql):
            deleted[row["id"]] = row[1:]

        if len(deleted) > 0:
            result["columns in deleted rows"] = row.keys()[1:]
            result["deleted"] = deleted

        # Detach the other database if needed
        if detach:
            self.system_db.detach(other_db)

        return result

    def _save(self):
        """Save our values in the database."""
        parameters = list(self._parameters)
        parameters.append(self.id)
        self.cursor.execute(
            "UPDATE cell SET a=?, b=?, c=?, alpha=?, beta=?, gamma=?" "WHERE id = ?",
            parameters,
        )
        self.db.commit()

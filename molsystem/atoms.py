# -*- coding: utf-8 -*-

"""A dictionary-like object for holding atoms

For efficiency the atom-data is stored as arrays -- numpy if possible, python
lists otherwise, along with metadata to define the attributes ("columns").

In some ways this a bit like pandas; however, we also need more control than a
simple pandas dataframe provides....
"""

import logging
from typing import Any, Dict, TypeVar

import numpy as np

import molsystem

System_tp = TypeVar("System_tp", "System", None)
Atoms_tp = TypeVar("Atoms_tp", "_Atoms", str, None)

logger = logging.getLogger(__name__)


class _Atoms(molsystem.table._Table):
    """The Atoms class holds arrays of attributes describing atoms

    Two attributes are required: the coordinates ('xyz') and atomic numbers
    ('atno'). Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Atoms can be added ('append') or removed ('delete').
    """

    def __init__(
        self,
        system: System_tp,
        table: str = 'atoms',
        other: Atoms_tp = None
    ) -> None:

        super().__init__(system, table, other)

        self._word = 'atoms'

        if self.table not in self.system:
            self._initialize()

    @property
    def n_atoms(self) -> int:
        """The number of atoms"""
        return self.n_rows

    @property
    def coordinate_type(self):
        """The type of coordinates: 'fractional' or 'Cartesian'"""
        return self._coordinate_type

    @coordinate_type.setter
    def coordinate_type(self, value):
        try:
            if 'fractional'.startswith(value.lower()):
                self._coordinate_type = 'fractional'
            elif 'cartesian'.startswith(value.lower()):
                self._coordinate_type = 'Cartesian'
            else:
                raise ValueError(
                    "The coordinate_type must be 'Cartesian' or 'fractional', "
                    "not '{}'".format(value)
                )
        except Exception:
            raise ValueError(
                "The coordinate_type must be 'Cartesian' or 'fractional', "
                "not '{}'".format(value)
            )

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more atoms

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Check keys and lengths of added atoms
        # if 'x' not in kwargs or 'y' not in kwargs or 'z' not in kwargs:
        #     raise KeyError("The coordinates are required!")

        # Need to handle the elements specially. Can give atomic numbers,
        # or symbols. By construction the references to elements are identical
        # to their atomic numbers.

        if 'symbol' in kwargs:
            symbols = kwargs.pop('symbol')
            kwargs['atno'] = self.to_atnos(symbols)

        super().append(**kwargs)

    def atoms(self, *args):
        """Return an iterator over the atoms."""
        return self.rows(*args)

    def to_atnos(self, symbols):
        """Convert element symbols to atomic numbers."""
        result = []

        parameters = [(x,) for x in symbols]
        for row in self.db.executemany(
            'SELECT atno WHERE symbol = ?', parameters
        ):
            result.append(row['atno'])

    def _initialize(self):
        """Set up the atom table. It needs to have a primary key
        to be used as a foreign key in the bonds table and elsewhere.
        """
        self.add_attribute('id', coltype='int', pk=True)
        self.add_attribute('atno', coltype='int', references='elements')
        self.add_attribute('x', coltype='float', index=True)
        self.add_attribute('y', coltype='float', index=True)
        self.add_attribute('z', coltype='float', index=True)


if __name__ == '__main__':  # pragma: no cover
    import numpy.random as nprand
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False) -> None:
        system = None
        atoms = _Atoms(system)
        x = nprand.uniform(low=0, high=100, size=nper)
        y = nprand.uniform(low=0, high=100, size=nper)
        z = nprand.uniform(low=0, high=100, size=nper)
        atno = np.array(nper * [6])
        with atoms as tmp:
            for i in range(0, nrepeat):
                tmp.append(x=x, y=y, z=z, atno=atno)
                x += 10.0
                y += 5.0
                z -= 3.0

    nrepeat = 1
    nper = 10
    t = timeit.timeit(
        "run(nper={}, nrepeat={})".format(nper, nrepeat),
        setup="from __main__ import run",
        timer=time.time,
        number=1
    )

    print("Creating {} atoms took {:.3f} s".format(nper * nrepeat, t))

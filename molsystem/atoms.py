# -*- coding: utf-8 -*-

import logging
import molsystem
import numpy as np


"""A dictionary-like object for holding atoms

For efficiency the atom-data is stored as arrays -- numpy if possible, python
lists otherwise, along with metadata to define the attributes ("columns").

In some ways this a bit like pandas; however, we also need more control than a
simple pandas dataframe provides....
"""

logger = logging.getLogger(__name__)

attributes = {
    'x': {
        'type': np.float64,
        'default': 0.0,
    },
    'y': {
        'type': np.float64,
        'default': 0.0,
    },
    'z': {
        'type': np.float64,
        'default': 0.0,
    },
    'atno': {
        'type': np.int8,
        'default': -1,
    },
    'name': {
        'type': np.object,
        'default': '',
    },
}


class Atoms(molsystem.Table):
    """The Atoms class holds arrays of attributes describing atoms

    Two attributes are required: the coordinates ('xyz') and atomic numbers
    ('atno'). Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Atoms can be added ('append') or removed ('delete').
    """

    def __init__(self, system=None):

        super().__init__(system)

        self._word = 'atoms'

        if not isinstance(system, molsystem.Atoms):
            self._private['system'] = system
            self._private['attributes'] = attributes
            for key in ('x', 'y', 'z', 'atno'):
                self.add_attribute(key)

    @property
    def n_atoms(self):
        """The number of atoms"""
        return self._private['n_rows']

    def append(self, **kwargs):
        """Append one or more atoms

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Check keys and lengths of added atoms
        if 'x' not in kwargs or 'y' not in kwargs or 'z' not in kwargs:
            raise KeyError("The coordinates are required!")

        super().append(**kwargs)


if __name__ == '__main__':  # pragma: no cover
    import numpy.random as nprand
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False):
        system = None
        atoms = Atoms(system)
        x = nprand.uniform(low=0, high=100, size=nper)
        y = nprand.uniform(low=0, high=100, size=nper)
        z = nprand.uniform(low=0, high=100, size=nper)
        atno = np.array(nper*[6])
        with atoms as tmp:
            if preallocate:
                tmp.allocate(nper*nrepeat)
            for i in range(0, nrepeat):
                tmp.append(x=x, y=y, z=z, atno=atno)
                x += 10.0
                y += 5.0
                z -= 3.0

    nrepeat = 1000
    nper = 100
    t = timeit.timeit("run(nper={}, nrepeat={})".format(nper, nrepeat),
                      setup="from __main__ import run",
                      timer=time.time, number=1)

    print("Creating {} atoms took {:.3f} s".format(nper*nrepeat, t))

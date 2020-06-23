# -*- coding: utf-8 -*-

import logging
from typing import Any, Dict, TypeVar

import numpy as np

import molsystem
"""A dictionary-like object for holding bonds

For efficiency the bond-data is stored as arrays -- numpy if possible, python
lists otherwise, along with metadata to define the attributes ("columns").

In some ways this a bit like pandas; however, we also need more control than a
simple pandas dataframe provides....
"""

System_tp = TypeVar("System_tp", "System", "Bonds", None)

logger = logging.getLogger(__name__)

attributes = {
    'i': {
        'type': np.int32,
        'default': -1,
    },
    'j': {
        'type': np.int32,
        'default': -1,
    },
    'order': {
        'type': np.int8,
        'default': 1,
    },
}


class Bonds(molsystem.Table):
    """The Bonds class holds arrays of attributes describing bonds

    Three attributes are required: the atoms 'i' and 'j' of the bond, and the
    bond order, which defaults to a single bond.

    Since there are two possible representations for a bond, i-j and j-i,
    the bonds are stored with i < j in order to make the indexing unique.

    Other attributes can be created, either from a list of predefined
    ones or by specifying the metadata required of an attribute. Attributes can
    also be removed. See the method 'add_attribute' for more detail.

    Bonds can be added ('append') or removed ('delete').
    """

    def __init__(self, system: System_tp = None) -> None:

        super().__init__(system)

        self._word = 'bonds'

        if not isinstance(system, molsystem.Bonds):
            self._private['system'] = system
            self._private['attributes'] = attributes
            for key in ('i', 'j', 'order'):
                self.add_attribute(key)

    @property
    def n_bonds(self) -> int:
        """The number of bonds"""
        return self._private['n_rows']

    def append(self, **kwargs: Dict[str, Any]) -> None:
        """Append one or more bonds

        The keys give the field for the data. If an existing field is not
        mentioned, then the default value is used, unless the default is None,
        in which case an error is thrown. It is an error if there is not a
        field corrresponding to a key.
        """

        # Check keys and lengths of added bonds
        if 'i' not in kwargs or 'j' not in kwargs:
            raise KeyError("The atoms i & j are required!")

        i = kwargs.pop('i')
        j = kwargs.pop('j')
        try:
            len_i = len(i)
        except TypeError:
            len_i = 1
            i = (i,)
        try:
            len_j = len(j)
        except TypeError:
            len_j = 1
            j = (j,)

        if len_i == 1:
            if len_j > 1:
                i = len_j * [i[0]]
        elif len_j == 1:
            j = len_i * [j[0]]
        elif len_i != len_j:
            raise IndexError(
                f'key "j" has the wrong number of values, {len_j}. '
                f'Should be 1 or the number of values in i, {len_i}.'
            )
        # Ensure that i < j

        if isinstance(i, np.ndarray):
            i = i.tolist()
        if isinstance(j, np.ndarray):
            j = j.tolist()

        i2 = []
        j2 = []
        for i_, j_ in zip(i, j):
            # will need to handle offsets here at some point
            if not isinstance(i_, int) or not isinstance(j_, int):
                print(i_)
                print(type(i_))
                raise TypeError(
                    "'i' and 'j', the atoms indices, must be integers"
                )
            if i_ < j_:
                i2.append(i_)
                j2.append(j_)
            else:
                i2.append(j_)
                j2.append(i_)

        super().append(i=i2, j=j2, **kwargs)


if __name__ == '__main__':  # pragma: no cover
    import timeit
    import time

    def run(nper=1000, nrepeat=100, preallocate=False) -> None:
        global i_atom, j_atom
        system = None
        bonds = Bonds(system)
        with bonds as tmp:
            if preallocate:
                tmp.allocate(n * nrepeat)
            for i in range(0, nrepeat):
                tmp.append(i=i_atom, j=j_atom)

    nrepeat = 1000
    nper = 100

    i_atom = []
    j_atom = []
    for i in range(0, nper, 3):
        if i > 2:
            i_atom.append(i - 3)
            j_atom.append(i)
        i_atom.append(i)
        j_atom.append(i + 1)
        i_atom.append(i)
        j_atom.append(i + 2)

    n = len(i_atom)

    t = timeit.timeit(
        "run(nrepeat={})".format(nrepeat),
        setup="from __main__ import run",
        timer=time.time,
        number=1
    )

    print("Creating {} bonds took {:.3f} s".format(n * nrepeat, t))

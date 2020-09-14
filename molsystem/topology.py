# -*- coding: utf-8 -*-

"""Topological methods for the system"""

import logging

logger = logging.getLogger(__name__)


class TopologyMixin:
    """A mixin for handling topology."""

    def find_molecules(self, configuration=None, as_indices=False):
        """Find the separate molecules in a system.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.

        Returns
        -------
        molecules : [][int]
            A list of lists of atom ids for the molecules
        """

        if configuration is None:
            configuration = self.current_configuration

        molecules = []

        atoms = self['atoms']
        atom_ids = atoms.atom_ids(configuration)
        n_atoms = len(atom_ids)

        if n_atoms == 0:
            return molecules

        to_index = {j: i for i, j in enumerate(atom_ids)}
        neighbors = self.bonded_neighbors(configuration)
        visited = [False] * n_atoms
        while True:
            # Find first atom not yet visited
            try:
                index = visited.index(False)
            except ValueError:
                break
            visited[index] = True
            i = atom_ids[index]
            atoms = [i]
            next_atoms = neighbors[i]
            while len(next_atoms) > 0:
                tmp = []
                for i in next_atoms:
                    if not visited[to_index[i]]:
                        atoms.append(i)
                        visited[to_index[i]] = True
                        tmp.extend(neighbors[i])
                next_atoms = tmp
            molecules.append(sorted(atoms))
        if as_indices:
            tmp = []
            for molecule in molecules:
                tmp.append([to_index[j] for j in molecule])
            return tmp
        else:
            return molecules

    def bonded_neighbors(self, configuration=None):
        """The atoms bonded to each atom in the system.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.

        Returns
        -------
        neighbors : {}[int]
            list of atom ids for each atom id
        """
        neighbors = {}

        atoms = self['atoms']
        bonds = self['bonds']
        n_atoms = atoms.n_atoms(configuration)

        if n_atoms == 0:
            return neighbors

        atom_ids = atoms.atom_ids(configuration)
        neighbors = {i: [] for i in atom_ids}

        if bonds.n_bonds(configuration) == 0:
            # No bonds, so just atoms....
            return neighbors

        for bond in bonds.bonds(configuration):
            i = bond['i']
            j = bond['j']
            neighbors[i].append(j)
            neighbors[j].append(i)

        return neighbors

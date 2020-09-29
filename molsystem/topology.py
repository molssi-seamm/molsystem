# -*- coding: utf-8 -*-

"""Topological methods for the system"""

import logging

from openbabel import openbabel

logger = logging.getLogger(__name__)


class TopologyMixin:
    """A mixin for handling topology."""

    def find_molecules(self, configuration=None, as_indices=False):
        """Find the separate molecules in a system.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.
        as_indices : bool = False
            Whether to return 0-based indices (True) or atom ids (False)

        Returns
        -------
        molecules : [[int]*n_molecules]
            A list of lists of atom ids or indices for the molecules
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
            return [[to_index[j] for j in js] for js in molecules]
        else:
            return molecules

    def bonded_neighbors(
        self, configuration=None, as_indices=False, first_index=0
    ):
        """The atoms bonded to each atom in the system.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.
        as_indices : bool = False
            Whether to return 0-based indices (True) or atom ids (False)
        first_index : int = 0
            The smallest index, e.g. 0 or 1

        Returns
        -------
        neighbors : {int: [int]} or [[int]] for indices
            list of atom ids for each atom id
        """
        neighbors = {}

        atoms = self['atoms']
        bonds = self['bonds']
        n_atoms = atoms.n_atoms(configuration)

        if n_atoms == 0:
            if as_indices:
                return []
            else:
                return neighbors

        atom_ids = atoms.atom_ids(configuration)
        neighbors = {i: [] for i in atom_ids}

        if bonds.n_bonds(configuration) > 0:
            for bond in bonds.bonds(configuration):
                i = bond['i']
                j = bond['j']
                neighbors[i].append(j)
                neighbors[j].append(i)

        if as_indices:
            # Convert to indices
            to_index = {j: i + first_index for i, j in enumerate(atom_ids)}
            result = [[]] * (n_atoms + first_index)
            for i, js in neighbors.items():
                result[to_index[i]] = sorted([to_index[j] for j in js])
            return result
        else:
            for i in neighbors:
                neighbors[i].sort()

            return neighbors

    def create_molecule_subsets(self, configuration=None):
        """Create a subset for each molecule in a configuration.

        By default they all reference an empty template 'all' of type
        'molecule'.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.

        Returns
        -------
        [int]
            The ids of the subsets, one per molecule.
        """

        if configuration is None:
            configuration = self.current_configuration

        # get the 'all/molecule' template
        templates = self['template']
        tid = templates.find('all', 'molecule', create=True)

        # Find the molecules and the create the subsets if they don't exist.
        molecules = self.find_molecules(configuration=configuration)

        # Remove any previous subsets for this configuration
        subsets = self['subset']
        sids = subsets.find(tid, configuration=configuration)
        if len(sids) > 0:
            subsets.delete(sids)

        # Now create the new set.
        sids = []
        for atom_ids in molecules:
            sid = subsets.create(
                tid, configuration=configuration, atoms=atom_ids
            )
            sids.append(sid)

        return sids

    def create_molecule_templates(
        self, configuration=None, create_subsets=True
    ):
        """Create a template for each unique molecule in a configuration.

        By default also create subsets linking each template to the atoms
        of the molecules in the system.

        Parameters
        ----------
        configuration : int = None
            The configuration to use, defaults to the current configuration.
        create_subsets : bool = True
            If true, create subsets linking the templates to the molecules.

        Returns
        -------
        [int] or [[int], [int]]
            The ids of the templates, or if create_subsets is True
            a two-element list containing the list of templates and
            list of subsets.
        """
        if configuration is None:
            configuration = self.current_configuration

        # Find the molecules
        molecules = self.find_molecules(configuration=configuration)
        n_molecules = len(molecules)

        # And the molecule each atom is in
        atom_to_molecule = {}
        for molecule, atoms in enumerate(molecules):
            for atom in atoms:
                atom_to_molecule[atom] = molecule

        # The bonds in each molecule
        bonds_per_molecule = [[] for i in range(n_molecules)]
        for bond in self.bonds.bonds(configuration=configuration):
            i = bond['i']
            j = bond['j']
            order = bond['bondorder']
            molecule = atom_to_molecule[i]
            bonds_per_molecule[molecule].append((i, j, order))

        # Get the canonical smiles for each molecule
        to_can = openbabel.OBConversion()
        to_can.SetOutFormat('can')
        to_smi = openbabel.OBConversion()
        to_smi.SetOutFormat('smi')
        ob_mol = openbabel.OBMol()
        ob_template = openbabel.OBMol()
        atnos = self.atoms.atomic_numbers(configuration)
        start = 0
        sids = {}
        tids = []
        for molecule, atoms in enumerate(molecules):
            to_index = {j: i for i, j in enumerate(atoms)}
            n_atoms = len(atoms)
            molecule_atnos = atnos[start:start + n_atoms]
            start += n_atoms

            ob_mol.Clear()
            for atom, atno in zip(atoms, molecule_atnos):
                ob_atom = ob_mol.NewAtom()
                ob_atom.SetAtomicNum(atno)
            bonds = []
            for i, j, order in bonds_per_molecule[molecule]:
                bonds.append((to_index[i], to_index[j], order))
                # 1-based indices in ob.
                ob_mol.AddBond(to_index[i] + 1, to_index[j] + 1, order)

            smiles = to_smi.WriteString(ob_mol).strip()
            canonical = to_can.WriteString(ob_mol).strip()

            # See if a molecule template with the canonical smiles exists
            if self.templates.exists(canonical, 'molecule'):
                tid = self.templates.find(canonical, 'molecule')
            else:
                tid = self.templates.create(
                    canonical, 'molecule', atnos=molecule_atnos, bonds=bonds
                )
            tids.append(tid)
            if create_subsets:
                tatom_ids = self.templateatoms.atom_ids(tid)
                if smiles != canonical:
                    # Need to reorder the atoms to match the template atoms

                    # Prepare the OB molecule for the template
                    ob_template.Clear()
                    for atno in self.templateatoms.atomic_numbers(tid):
                        ob_atom = ob_template.NewAtom()
                        ob_atom.SetAtomicNum(atno)

                    tatom_ids = self.templateatoms.atom_ids(tid)
                    to_index = {j: i for i, j in enumerate(tatom_ids)}
                    for row in self.templatebonds.bonds(tid):
                        i = to_index[row['i']]
                        j = to_index[row['j']]
                        order = row['bondorder']
                        ob_template.AddBond(i + 1, j + 1, order)

                    # Get the mapping from template to molecule
                    query = openbabel.CompileMoleculeQuery(ob_template)
                    mapper = openbabel.OBIsomorphismMapper.GetInstance(query)
                    mapping = openbabel.vpairUIntUInt()
                    mapper.MapFirst(ob_mol, mapping)
                    tmp = [atoms[j] for i, j in mapping]
                    atoms = tmp
                sid = self.subsets.create(tid, configuration, atoms, tatom_ids)
                if tid not in sids:
                    sids[tid] = []
                sids[tid].append(sid)

        if create_subsets:
            return tids, sids
        else:
            return tids

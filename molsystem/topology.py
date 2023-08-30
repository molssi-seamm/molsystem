# -*- coding: utf-8 -*-

"""Topological methods for the system"""

import logging
from math import floor
import pprint  # noqa: F401

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise

logger = logging.getLogger(__name__)


class TopologyMixin:
    """A mixin for handling topology in a configuration."""

    def find_molecules(self, as_indices=False):
        """Find the separate molecules.

        Parameters
        ----------
        as_indices : bool = False
            Whether to return 0-based indices (True) or atom ids (False)

        Returns
        -------
        molecules : [[int]*n_molecules]
            A list of lists of atom ids or indices for the molecules
        """

        molecules = []

        n_atoms = self.atoms.n_atoms

        if n_atoms == 0:
            return molecules

        if not as_indices and self.symmetry.n_symops > 1:
            raise RuntimeError(
                "Cannot return atom ids for bonded_neighbors when there is symmetry"
            )
        neighbors = self.bonded_neighbors(as_indices=True)
        visited = [False] * n_atoms
        while True:
            # Find first atom not yet visited
            try:
                i = visited.index(False)
            except ValueError:
                break
            visited[i] = True
            atoms = [i]
            next_atoms = neighbors[i]
            while len(next_atoms) > 0:
                tmp = []
                for i in next_atoms:
                    if not visited[i]:
                        atoms.append(i)
                        visited[i] = True
                        tmp.extend(neighbors[i])
                next_atoms = tmp
            molecules.append(sorted(atoms))
        if as_indices:
            return molecules
        else:
            atom_ids = self.atoms.ids
            to_id = {i: j for i, j in enumerate(atom_ids)}
            return [[to_id[j] for j in js] for js in molecules]

    def bonded_neighbors(self, as_indices=False, first_index=0):
        """The atoms bonded to each atom in the system.

        Parameters
        ----------
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

        symmetry = self.symmetry
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = atoms.n_atoms

        if n_atoms == 0:
            if as_indices:
                return []
            else:
                return neighbors

        if self.symmetry.n_symops > 1:
            # Symmetry involved...
            if not as_indices:
                raise RuntimeError(
                    "Cannot return atom ids for bonded_neighbors when there is symmetry"
                )
            pairs = symmetry.bond_atoms
            neighbors = [[] for i in range(n_atoms)]
            for i, j in pairs:
                neighbors[i].append(j)
                neighbors[j].append(i)
            for js in neighbors:
                js.sort()
            return neighbors
        else:
            atom_ids = atoms.ids
            neighbors = {i: [] for i in atom_ids}

            if bonds.n_bonds > 0:
                for bond in bonds.bonds():
                    i = bond["i"]
                    j = bond["j"]
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

    def create_molecule_subsets(self):
        """Create a subset for each molecule in a configuration.

        Returns
        -------
        [int]
            The ids of the subsets, one per molecule.
        """
        # Find the molecules and the create the subsets if they don't exist.
        molecules = self.find_molecules()

        # Remove any previous subsets for this configuration
        subsets = self["subset"]
        tid = 1
        sids = subsets.find(tid)
        if len(sids) > 0:
            subsets.delete(sids)

        # Now create the new set.
        sids = []
        for atom_ids in molecules:
            sid = subsets.create(tid, atoms=atom_ids)
            sids.append(sid)

        return sids

    def create_molecule_templates(self, full_templates=True, create_subsets=True):
        """Create a template for each unique molecule in a configuration.

        By default also create subsets linking each template to the atoms
        of the molecules in the system.

        Parameters
        ----------
        full_templates : bool = True
            If true, create full templates by creating systems for the
            molecules.
        create_subsets : bool = True
            If true, create subsets linking the templates to the molecules.

        Returns
        -------
        [int] or [[int], [int]]
            The ids of the templates, or if create_subsets is True
            a two-element list containing the list of templates and
            list of subsets.
        """
        templates = self.system_db.templates

        # Find the molecules
        molecules = self.find_molecules()
        n_molecules = len(molecules)

        # And the molecule each atom is in
        atom_to_molecule = {}
        for molecule, atoms in enumerate(molecules):
            for atom in atoms:
                atom_to_molecule[atom] = molecule

        # The bonds in each molecule
        bonds_per_molecule = [[] for i in range(n_molecules)]
        for bond in self.bonds.bonds():
            i = bond["i"]
            j = bond["j"]
            order = bond["bondorder"]
            molecule = atom_to_molecule[i]
            bonds_per_molecule[molecule].append((i, j, order))

        # Get the canonical smiles for each molecule
        to_can = openbabel.OBConversion()
        to_can.SetOutFormat("can")
        ob_mol = openbabel.OBMol()
        ob_template = openbabel.OBMol()
        atnos = self.atoms.atomic_numbers
        xyzs = self.atoms.get_coordinates(fractionals=False)

        atom_index = {j: i for i, j in enumerate(self.atoms.ids)}

        new_subsets = {}
        sids = {}
        new_templates = []
        tids = []
        for molecule, atoms in enumerate(molecules):
            to_index = {j: i for i, j in enumerate(atoms)}
            molecule_atnos = [atnos[atom_index[i]] for i in atoms]

            ob_mol.Clear()
            for atom, atno in zip(atoms, molecule_atnos):
                ob_atom = ob_mol.NewAtom()
                ob_atom.SetAtomicNum(atno)
            bonds = []
            for i, j, order in bonds_per_molecule[molecule]:
                bonds.append((to_index[i], to_index[j], order))
                # 1-based indices in ob.
                ob_mol.AddBond(to_index[i] + 1, to_index[j] + 1, order)

            canonical = to_can.WriteString(ob_mol).strip()

            if full_templates:
                # See if a molecule template with the canonical smiles exists
                if templates.exists(canonical, "molecule"):
                    template = templates.get(canonical, category="molecule")
                else:
                    # Create a new system & configuration for the template
                    system_name = "template system " + canonical
                    if not self.system_db.system_exists(system_name):
                        system = self.system_db.create_system(
                            system_name, make_current=False
                        )
                        configuration = system.create_configuration(
                            canonical, make_current=False
                        )
                        cid = configuration.id

                        kwargs = {}
                        kwargs["atno"] = molecule_atnos
                        molecule_xyzs = [xyzs[atom_index[i]] for i in atoms]
                        kwargs["x"] = [x for x, y, z in molecule_xyzs]
                        kwargs["y"] = [y for x, y, z in molecule_xyzs]
                        kwargs["z"] = [z for x, y, z in molecule_xyzs]

                        ids = configuration.atoms.append(**kwargs)

                        kwargs = {}
                        kwargs["i"] = [ids[x] for x, _, _ in bonds]
                        kwargs["j"] = [ids[x] for _, x, _ in bonds]
                        kwargs["bondorder"] = [x for _, _, x in bonds]

                        configuration.bonds.append(**kwargs)

                        template = templates.create(
                            canonical, category="molecule", configuration=cid
                        )
            else:
                if templates.exists(canonical, "molecule"):
                    template = templates.get(canonical, category="molecule")
                else:
                    template = templates.create(canonical, category="molecule")

            if template.id not in tids:
                tids.append(template.id)
                new_templates.append(template)

            if create_subsets:
                if full_templates:
                    # Need to reorder the atoms to match the template atoms
                    # Prepare the OB molecule for the template
                    ob_template.Clear()
                    for atno in template.atoms.atomic_numbers:
                        ob_atom = ob_template.NewAtom()
                        ob_atom.SetAtomicNum(atno)

                    tatom_ids = template.atoms.ids
                    to_index = {j: i for i, j in enumerate(tatom_ids)}
                    for row in template.bonds.bonds():
                        i = to_index[row["i"]]
                        j = to_index[row["j"]]
                        order = row["bondorder"]
                        ob_template.AddBond(i + 1, j + 1, order)

                    # Get the mapping from template to molecule
                    query = openbabel.CompileMoleculeQuery(ob_template)
                    mapper = openbabel.OBIsomorphismMapper.GetInstance(query)
                    mapping = openbabel.vpairUIntUInt()
                    mapper.MapFirst(ob_mol, mapping)
                    reordered_atoms = [atoms[j] for i, j in mapping]

                    if len(reordered_atoms) == 0:
                        logger.warning(
                            "There are no atoms in the subset for this molecule!"
                        )
                    else:
                        subset = self.subsets.create(
                            template, reordered_atoms, tatom_ids
                        )
                else:
                    subset = self.subsets.create(template, atoms)

                tid = template.id
                if tid not in sids:
                    sids[tid] = []
                    new_subsets[tid] = []
                sids[tid].append(subset.id)
                new_subsets[tid].append(subset)

        if create_subsets:
            return new_templates, new_subsets
        else:
            return new_templates

    def get_molecule_smiles(self):
        """Return the a list of the canonical SMILES for each molecule..

        Returns
        -------
        [str]
            The canonical SMILES for each molecule, in order that they are found.
        """
        # Find the molecules
        molecules = self.find_molecules()
        n_molecules = len(molecules)

        # And the molecule each atom is in
        atom_to_molecule = {}
        for molecule, atoms in enumerate(molecules):
            for atom in atoms:
                atom_to_molecule[atom] = molecule

        # The bonds in each molecule
        bonds_per_molecule = [[] for i in range(n_molecules)]
        for bond in self.bonds.bonds():
            i = bond["i"]
            j = bond["j"]
            order = bond["bondorder"]
            molecule = atom_to_molecule[i]
            bonds_per_molecule[molecule].append((i, j, order))

        # Get the canonical smiles for each molecule
        to_can = openbabel.OBConversion()
        to_can.SetOutFormat("can")
        ob_mol = openbabel.OBMol()
        atnos = self.atoms.atomic_numbers

        atom_index = {j: i for i, j in enumerate(self.atoms.ids)}

        SMILES = []
        for molecule, atoms in enumerate(molecules):
            to_index = {j: i for i, j in enumerate(atoms)}
            molecule_atnos = [atnos[atom_index[i]] for i in atoms]

            ob_mol.Clear()
            for atom, atno in zip(atoms, molecule_atnos):
                ob_atom = ob_mol.NewAtom()
                ob_atom.SetAtomicNum(atno)
            bonds = []
            for i, j, order in bonds_per_molecule[molecule]:
                bonds.append((to_index[i], to_index[j], order))
                # 1-based indices in ob.
                ob_mol.AddBond(to_index[i] + 1, to_index[j] + 1, order)

            canonical = to_can.WriteString(ob_mol).strip()
            SMILES.append(canonical)
        return SMILES

    def reimage_bonded_atoms(self, reimage_molecules=True):
        """Ensure that the atoms in a molecule are "near" each other.

        In a periodic system atoms can be translated by a unit cell without changing
        anything. However the bonds in a molecule need to account for such translations
        of atoms otherwise a bond may appear to be very long.

        It is often convenient physically move atoms by the correct cell translations
        to bring the atoms of a molecule close to each other. That is what this method
        does, moving all the atoms close to the first. Optionally the molecule is
        also moved to bring its geometric center into the primary unit cell, i.e. with
        fractional coordinates in the range [0..1).

        Parameters
        ----------
        reimage_molecules : bool = True
            Whether to move molecules into the primary unit cell.

        Returns
        -------
        bool
            True if the coordinates were changed.
        """
        if self.periodicity == 0:
            # Nothing to do.
            return False

        logger.log(0, f"Configuration {self.id} {self}")
        logger.log(0, self.atoms)
        logger.log(0, 10 * "+")
        # The bonds from each atom
        bonds = self.bonded_neighbors(as_indices=True)

        logger.debug("Bonds:")
        logger.debug(pprint.pformat(bonds))
        logger.debug("")

        # And coordinates as fractionals
        tmp = self.atoms.get_coordinates(fractionals=True)
        xyz = [[x, y, z] for x, y, z in tmp]

        visited = [False] * len(xyz)

        # Find the molecules
        molecules = self.find_molecules(as_indices=True)

        # Work through the atoms, see if bonded atoms need to be moved
        changed = False
        for atoms in molecules:
            logger.debug("")
            logger.debug(f"{atoms=}")
            # if len(atoms) > 1:
            #     atom = sorted(atoms)[0]
            #     logger.debug(f" --> {atom}")
            #     if self._image_helper(bonds, xyz, atom):
            #         changed = True

            changed = self._image_helper2(bonds, xyz, visited, atoms[0], changed)
        # Check temporarily.
        for iatom, jatoms in enumerate(bonds):
            xyzi = xyz[iatom]
            for jatom in jatoms:
                xyzj = xyz[jatom]
                for vi, vj in zip(xyzi, xyzj):
                    delta = floor(vj - vi + 0.5)
                    if delta != 0:
                        print(
                            f"Warning: {iatom} - {jatom} delta = {vj:.3f} - {vi:.3f} = "
                            f"{delta}"
                        )

        # Move the center of molecules into the primary cell, if requested.
        if reimage_molecules:
            for atoms in molecules:
                cx = cy = cz = 0.0
                for atom in atoms:
                    x, y, z = xyz[atom]
                    cx += x
                    cy += y
                    cz += z
                n = len(atoms)
                dx = int(cx / n)
                dy = int(cy / n)
                dz = int(cz / n)
                if dx != 0 or dy != 0 or dz != 0:
                    changed = True
                    for atom in atoms:
                        x, y, z = xyz[atom]
                        xyz[atom] = [x - dx, y - dy, z - dz]
        if changed:
            self.atoms.set_coordinates(xyz, fractionals=True)

        return changed

    def _image_helper(self, bonds, xyz, atom):
        """Translate all the bonded neighbors with larger index.

        Parameters
        ----------
        bonds : [[int]]
            The bonded neighbors for each atom.
        xyz : [[float * 3]]
            The fractional coordinates of the atoms.
        """
        logger.debug(f"       helper {atom}: {bonds[atom]}")
        changed = False
        xyzi = xyz[atom]
        for jatom in bonds[atom]:
            logger.debug(f"\t\t{atom=} {jatom}")
            if jatom < atom:
                continue
            xyzj = xyz[jatom]
            new_xyz = []
            shift = False
            for vi, vj in zip(xyzi, xyzj):
                delta = floor(vj - vi + 0.5)
                if delta == 0:
                    new_xyz.append(vj)
                else:
                    new_xyz.append(vj - delta)
                    shift = True
                logger.log(
                    0, f"{atom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
            if shift:
                changed = True
                logger.debug(
                    f"{atom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
                xyz[jatom] = new_xyz

        for jatom in bonds[atom]:
            if jatom < atom:
                continue
            if self._image_helper(bonds, xyz, jatom):
                changed = True

        return changed

    def reimage_molecules(self):
        """Reimage molecules into the primary unit cell.

        The molecules are moved to bring their geometric center into the primary unit
        cell, i.e. with fractional coordinates in the range [0..1).

        Returns
        -------
        bool
            True if the coordinates were changed.
        """
        if self.periodicity == 0:
            # Nothing to do.
            return False

        # And coordinates as fractionals
        xyz = self.atoms.get_coordinates(fractionals=True)

        # Find the molecules
        molecules = self.find_molecules(as_indices=True)

        # Move the centor of molecules into the primary cell, if requested.
        changed = False
        for atoms in molecules:
            cx = cy = cz = 0.0
            for atom in atoms:
                x, y, z = xyz[atom]
                cx += x
                cy += y
                cz += z
            n = len(atoms)
            dx = int(cx / n)
            dy = int(cy / n)
            dz = int(cz / n)
            if dx != 0 or dy != 0 or dz != 0:
                changed = True
                for atom in atoms:
                    x, y, z = xyz[atom]
                    xyz[atom] = [x - dx, y - dy, z - dz]
        if changed:
            self.atoms.set_coordinates(xyz, fractionals=True)

        return changed

    def _image_helper2(self, bonds, xyz, visited, iatom, changed):
        visited[iatom] = True
        xyzi = xyz[iatom]
        for jatom in bonds[iatom]:
            logger.debug(f"\t\t{iatom=} {jatom}")
            if visited[jatom]:
                continue
            new_xyz = []
            xyzj = xyz[jatom]
            shift = False
            for vi, vj in zip(xyzi, xyzj):
                delta = floor(vj - vi + 0.5)
                if delta == 0:
                    new_xyz.append(vj)
                else:
                    new_xyz.append(vj - delta)
                    changed = True
                    shift = True
            if shift:
                logger.debug(
                    f"{iatom:3d} - {jatom:3d}: {xyzi} {xyzj} --> {new_xyz} {shift}"
                )
            xyz[jatom] = new_xyz
            changed = self._image_helper2(bonds, xyz, visited, jatom, changed)
        return changed

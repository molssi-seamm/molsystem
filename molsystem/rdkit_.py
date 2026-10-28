# -*- coding: utf-8 -*-

"""Interface to RDKit."""

import logging
import json
import pprint
import string

try:
    from rdkit import Chem
except ModuleNotFoundError:
    print(
        "Please install rdkit using conda:\n" "     conda install -c conda-forge rdkit"
    )
    raise

from seamm_util import CompactJSONEncoder

logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)


# Pauling electronegativities, to place the charges that the net charge needs:
# negative charges on the more, positive ones on the less electronegative atoms.
_electronegativity = {
    1: 2.20,
    5: 2.04,
    6: 2.55,
    7: 3.04,
    8: 3.44,
    9: 3.98,
    14: 1.90,
    15: 2.19,
    16: 2.58,
    17: 3.16,
    33: 2.18,
    34: 2.55,
    35: 2.96,
    53: 2.66,
}

# Give up on the search for charges and radicals beyond this many partial
# assignments, e.g. a large system whose bonds were perceived without orders.
_max_states = 20000


def _charge_cost(atno, q, metal=False):
    """The cost of putting the formal charge ``q`` on an atom.

    Charges are costly, negative ones less so on electronegative atoms and
    positive ones on electropositive atoms, and charges on carbon more, since
    carbocations and carbanions are rarer than ions of heteroatoms. A metal is
    better as a cation than neutral, e.g. Li+ in an electrolyte.
    """
    if metal:
        if q > 0:
            return 0.1 * q
        return 1.0 + 2.0 * abs(q)
    if q == 0:
        return 0.0
    en = _electronegativity.get(atno, 2.5)
    if q < 0:
        cost = abs(q) * (1.0 + 0.1 * (4.0 - en))
    else:
        cost = q * (1.0 + 0.1 * en)
    if atno == 6:
        cost += abs(q)
    return cost


def _candidates(atno, v_sigma, aromatic, charge=None):
    """The (formal charge, radical electrons, pi) an atom can have.

    The valence an atom with formal charge q may have is that of the
    isoelectronic element, as RDKit has it: N+ like C, O- like F, B- like C,
    P- like S. An atom short of its lowest valence has the difference as
    radical electrons. Atoms in aromatic rings may take part in a ring double
    bond (pi = 1), like the carbons of benzene, or not (pi = 0), like the N of
    pyrrole or a carbanion in a ring.

    Parameters
    ----------
    atno : int
        The atomic number.
    v_sigma : int
        The valence from the bonds, counting each aromatic bond as single.
    aromatic : bool
        Whether the atom has aromatic bonds.
    charge : int = None
        The formal charge, if known; otherwise charges from -2 to +2 are tried.

    Returns
    -------
    [(int, int, int)]
        The possible (charge, radicals, pi).
    """
    table = Chem.GetPeriodicTable()
    result = []
    for q in (0, -1, 1, -2, 2) if charge is None else (charge,):
        z = atno - q
        if z < 1 or z > 118:
            continue
        valences = list(table.GetValenceList(z))
        if -1 in valences:
            # Metals have any valence and no radicals in RDKit
            if q == 0:
                result.append((0, 0, 0))
            continue
        for pi in (0, 1) if aromatic else (0,):
            v = v_sigma + pi
            allowed = [x for x in valences if x >= v]
            if len(allowed) == 0:
                continue
            radicals = allowed[0] - v
            if radicals <= 4:
                result.append((q, radicals, pi))
    return result


def _charges_and_radicals(
    atnos, bonds, charge=None, n_unpaired=None, fixed=None, prune=True
):
    """The formal charges and radical electrons of the atoms of a structure.

    All the hydrogens are explicit, so an atom short of its valence is a radical,
    an anion or a cation, depending on the net charge and the number of unpaired
    electrons. The choice is the one that matches the net charge and the number of
    unpaired electrons, with the fewest charges, negative charges on electronegative
    atoms and positive ones on electropositive atoms; e.g. the ethyl radical
    [CH2]C, the ethyl cation [CH2+]C, nitromethane C[N+](=O)[O-], sulfuric acid
    O=S(=O)(O)O and PF6- with the charge on P. The ring double bonds of an aromatic
    system pair its atoms, so an even number of them must take part.

    Parameters
    ----------
    atnos : [int]
        The atomic numbers.
    bonds : [(int, int, int)]
        The bonds as (i, j, bond order), with indices into atnos and 5 for an
        aromatic bond.
    charge : int = None
        The net charge, or None if not known.
    n_unpaired : int = None
        The number of unpaired electrons, or None if not known.
    fixed : [int] = None
        Known formal charges of the atoms.
    prune : bool = True
        Keep neutral, closed-shell atoms outside aromatic rings so, which keeps the
        search small for large systems but misses e.g. radical cations like C[NH2+].

    Returns
    -------
    ([int], [int])
        The formal charges and the numbers of radical electrons, or (None, None)
        if there are too many possibilities to search.
    """
    n = len(atnos)
    v_sigma = [0] * n
    aromatic_neighbors = [[] for _ in range(n)]
    for i, j, order in bonds:
        if order == 5:
            v_sigma[i] += 1
            v_sigma[j] += 1
            aromatic_neighbors[i].append(j)
            aromatic_neighbors[j].append(i)
        else:
            v_sigma[i] += order
            v_sigma[j] += order

    # The aromatic systems, each a set of atoms joined by aromatic bonds
    system = [None] * n
    n_systems = 0
    for start in range(n):
        if system[start] is None and len(aromatic_neighbors[start]) > 0:
            stack = [start]
            system[start] = n_systems
            while stack:
                i = stack.pop()
                for j in aromatic_neighbors[i]:
                    if system[j] is None:
                        system[j] = n_systems
                        stack.append(j)
            n_systems += 1

    # Metals, which RDKit allows any valence, may be ions, e.g. Li+
    table = Chem.GetPeriodicTable()
    metals = {z for z in set(atnos) if -1 in list(table.GetValenceList(z))}

    candidates = []
    for i, atno in enumerate(atnos):
        tmp = _candidates(
            atno,
            v_sigma[i],
            system[i] is not None,
            charge=None if fixed is None else fixed[i],
        )
        if len(tmp) == 0:
            # A valence nothing fits: leave the atom as it is, for RDKit to judge
            tmp = [(0 if fixed is None else fixed[i], 0, 0)]
        elif prune and system[i] is None and (0, 0, 0) in tmp and atno not in metals:
            # A neutral, closed-shell atom outside aromatic rings stays so, which
            # keeps the search small for large systems, e.g. a box of water
            tmp = [(0, 0, 0)]
        candidates.append(tmp)

    radical_cost = 0.0 if n_unpaired is not None else 1.5

    # Dynamic programming over the atoms, the aromatic systems last and each
    # together, so that the parity of its ring double bonds can be checked at its
    # end. A state is (sum of charges, sum of radicals, parity) -> (cost, choices),
    # with the choices as a linked list.
    order = sorted(range(n), key=lambda i: -1 if system[i] is None else system[i])
    states = {(0, 0, 0): (0.0, None)}
    for k, i in enumerate(order):
        new_states = {}
        for (sq, sr, parity), (cost, choices) in states.items():
            for q, r, pi in candidates[i]:
                key = (sq + q, sr + r, parity ^ pi)
                value = (
                    cost
                    + _charge_cost(atnos[i], q, metal=atnos[i] in metals)
                    + radical_cost * r
                )
                if key not in new_states or value < new_states[key][0]:
                    new_states[key] = (value, ((i, q, r), choices))
        states = new_states
        last = k + 1 == n or system[order[k + 1]] != system[i]
        if system[i] is not None and last:
            # The end of an aromatic system: its ring double bonds must pair up
            states = {
                (sq, sr, 0): value
                for (sq, sr, parity), value in states.items()
                if parity == 0
            }
            if len(states) == 0:
                # No way to pair them; carry on without the check
                states = new_states
                states = {(sq, sr, 0): value for (sq, sr, _), value in states.items()}
        if len(states) > _max_states:
            return None, None

    def badness(item):
        (sq, sr, _), (cost, _) = item
        return (
            0 if charge is None else abs(sq - charge),
            0 if n_unpaired is None else abs(sr - n_unpaired),
            cost,
        )

    _, (_, choices) = min(states.items(), key=badness)
    charges = [0] * n
    radicals = [0] * n
    while choices is not None:
        (i, q, r), choices = choices
        charges[i] = q
        radicals[i] = r
    return charges, radicals


references = [
    """\
@misc{rdkit,
    title = {RDKit: Open-source cheminformatics.},
    month = {$month},
    year = {$year},
    author = {Greg Landrum, et al.},
    url = {https://www.rdkit.org},
    doi = {10.5281/zenodo.591637},
    version = {$version}
}
    """,
]


def rdkit_version():
    """Return the RDKit version."""
    return Chem.rdBase.rdkitVersion


def rdkit_citations():
    """Return the BibTeX citations for RDKit"""
    citations = []
    try:
        version = rdkit_version()
        tmp = version.split(".")
        year = tmp[0]
        if len(tmp) > 1:
            month = tmp[1]
        else:
            month = "?"

        for reference in references:
            template = string.Template(reference)
            citation = template.substitute(month=month, version=version, year=year)
            citations.append(citation)
    except Exception:
        pass
    return citations


class RDKitMixin:
    """A mixin for handling RDKit via its Python interface."""

    def coordinates_from_RDKMol(self, mol_or_conformer):
        """Update the coordinates from an RDKMol.

        Parameters
        ----------
        mol_or_conformer : rdkit.Chem.Mol or rdkit.Chem.rdchem.Conformer
            The RDKit molecule or conformer object
        """
        if isinstance(mol_or_conformer, Chem.Mol):
            if mol_or_conformer.GetNumConformers() == 0:
                raise ValueError("RDKit molecule has no conformers")
            conformer = mol_or_conformer.GetConformer()
        elif isinstance(mol_or_conformer, Chem.rdchem.Conformer):
            conformer = mol_or_conformer
        else:
            raise TypeError(
                f"Don't recognize molecule/conformer argument {type(mol_or_conformer)}"
            )
        self.atoms.coordinates = conformer.GetPositions()

    def coordinates_to_RDKMol(self, rdkmol):
        """Update the coordinates of an RDKMol from this configuration.

        Parameters
        ----------
        rdkmol : rdkit.Chem.RWMol
            The RDKMol to use.
        """
        rdkconf = rdkmol.GetConformers()[0]
        rdkconf.SetPositions(self.atoms.coordinates)

    def to_RDKMol(self, properties=None):
        """Return an RDKMol object for the configuration, template, or subset."""
        index = {}
        indices = []
        rdk_mol = Chem.RWMol()
        for atno, _id in zip(self.atoms.atomic_numbers, self.atoms.ids):
            atom = Chem.Atom(atno)
            idx = rdk_mol.AddAtom(atom)
            index[_id] = idx
            indices.append(idx)

        bond_types = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            5: Chem.BondType.AROMATIC,
        }
        for row in self.bonds.bonds():
            if row["bondorder"] > 0:
                rdk_mol.AddBond(
                    index[row["i"]], index[row["j"]], bond_types[row["bondorder"]]
                )

        # The formal charges and radicals. The hydrogens are all explicit, so RDKit
        # must not add any: an atom short of its valence is a radical or an ion.
        atnos = self.atoms.atomic_numbers
        bonds = [
            (index[row["i"]], index[row["j"]], row["bondorder"])
            for row in self.bonds.bonds()
            if row["bondorder"] in bond_types
        ]
        if self.__class__.__name__ == "_Configuration":
            charge = self.charge
            n_unpaired = self.spin_multiplicity - 1
        else:
            charge = None
            n_unpaired = None
        fixed = None
        if "formal_charge" in self.atoms:
            fixed = self.atoms.get_column_data("formal_charge")
            if charge is not None and sum(fixed) != charge:
                fixed = None
        # The formal charges stored with the atoms are kept unless others fit the
        # net charge and spin better, e.g. after the bonds have changed.
        best = None
        tries = [(None, True), (None, False)]
        if fixed is not None:
            tries.insert(0, (fixed, True))
        for known, prune in tries:
            charges, radicals = _charges_and_radicals(
                atnos,
                bonds,
                charge=charge,
                n_unpaired=n_unpaired,
                fixed=known,
                prune=prune,
            )
            if charges is None:
                continue
            miss = (
                0 if charge is None else abs(sum(charges) - charge),
                0 if n_unpaired is None else abs(sum(radicals) - n_unpaired),
            )
            if best is None or miss < best[0]:
                best = (miss, charges, radicals)
            if miss == (0, 0):
                break
        if best is not None:
            _, charges, radicals = best
        else:
            charges, radicals = fixed, None
        for idx, at in enumerate(rdk_mol.GetAtoms()):
            at.SetNoImplicit(True)
            if charges is not None:
                at.SetFormalCharge(charges[idx])
            if radicals is not None:
                at.SetNumRadicalElectrons(radicals[idx])
        for bond in rdk_mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetIsAromatic(True)
                bond.GetBeginAtom().SetIsAromatic(True)
                bond.GetEndAtom().SetIsAromatic(True)

        natom = self.atoms.n_atoms
        conf = Chem.Conformer(natom)
        for idx, xyz in zip(indices, self.atoms.coordinates):
            conf.SetAtomPosition(idx, xyz)

        rdk_mol.AddConformer(conf)
        try:
            Chem.SanitizeMol(rdk_mol)
        except Chem.KekulizeException as e:
            logger.warning(f"Kekulization failed: {e}")
            flags = (
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            Chem.rdmolops.SanitizeMol(rdk_mol, sanitizeOps=flags)
        except Chem.AtomValenceException as e:
            # e.g. CH4+, which RDKit cannot write
            logger.warning(f"RDKit does not accept the valences: {e}")
            flags = (
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            )
            Chem.rdmolops.SanitizeMol(rdk_mol, sanitizeOps=flags)

        # Add the net charge and spin multiplicity as properties for configurations
        if self.__class__.__name__ == "_Configuration":
            rdk_mol.SetIntProp("SEAMM|net charge|int|", self.charge)
            rdk_mol.SetIntProp("SEAMM|spin multiplicity|int|", self.spin_multiplicity)
            rdk_mol.SetProp(
                "SEAMM|XYZ|json|",
                json.dumps(self.coordinates, indent=4, cls=CompactJSONEncoder),
            )
            if self.system.name is not None:
                rdk_mol.SetProp("SEAMM|system name|str|", self.system.name)
            if self.name is not None:
                rdk_mol.SetProp("SEAMM|configuration name|str|", self.name)

        if properties is not None:
            data = self.properties.get(properties, include_system_properties=True)
            for _property, value in data.items():
                _type = self.properties.type(_property)
                units = self.properties.units(_property)
                value = value["value"]
                key = f"SEAMM|{_property}|{_type}|"
                if units is not None and units != "":
                    key += units

                if _type == "int":
                    rdk_mol.SetIntProp(key, value)
                elif _type == "float":
                    rdk_mol.SetDoubleProp(key, value)
                elif _type == "json":
                    rdk_mol.SetProp(key, json.dumps(value))
                else:
                    rdk_mol.SetProp(key, str(value))

        return rdk_mol

    def from_RDKMol(
        self,
        mol_or_conformer,
        properties="all",
        atoms=True,
        coordinates=True,
        bonds=True,
    ):
        """Transform an RDKit molecule into the current object.

        Parameters
        ----------
        mol_or_conformer : rdkit.Chem.Mol or rdkit.Chem.rdchem.Conformer
            The RDKit molecule or conformer object

        properties : str = "all"
            Whether to include all properties or none

        atoms : bool = True
            Recreate the atoms

        coordinates : bool = True
            Update the coordinates

        bonds : bool = True
            Recreate the bonds from the RDKit molecule

        Returns
        -------
        (str, str) : system_name, configuration_name
        """
        system_name = None
        configuration_name = None

        if isinstance(mol_or_conformer, Chem.Mol):
            mol = mol_or_conformer
            if mol.GetNumConformers() == 0:
                raise ValueError("RDKit molecule has no conformers")
            conformer = mol.GetConformer()
        elif isinstance(mol_or_conformer, Chem.rdchem.Conformer):
            conformer = mol_or_conformer
            mol = conformer.GetOwningMol()
        else:
            raise RuntimeError(
                f"Don't recognize molecule/conformer argument {type(mol_or_conformer)}"
            )

        atnos = []
        qs = []
        for rdk_atom in mol.GetAtoms():
            atnos.append(rdk_atom.GetAtomicNum())
            qs.append(rdk_atom.GetFormalCharge())
            logger.debug(f"atom {atnos}")

        Xs = []
        Ys = []
        Zs = []
        for atom_coords in conformer.GetPositions():
            Xs.append(atom_coords[0])
            Ys.append(atom_coords[1])
            Zs.append(atom_coords[2])
            logger.debug(f"{atom_coords[0]} {atom_coords[1]} {atom_coords[2]}")

        Is = []
        Js = []
        BondOrders = []
        bond_types = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 5,
        }
        for rdk_bond in mol.GetBonds():
            i = rdk_bond.GetBeginAtom().GetIdx()
            j = rdk_bond.GetEndAtom().GetIdx()
            bondorder = bond_types[rdk_bond.GetBondType()]
            Is.append(i)
            Js.append(j)
            BondOrders.append(bondorder)
            logger.debug(f"bond {i} - {j} {bondorder}")

        if atoms:
            self.clear()
            try:
                # The formal charges, as from_OBMol keeps them: analyses of an MD
                # without fixed charges (an MLFF) find the ions from them
                if any(q != 0 for q in qs):
                    if "formal_charge" not in self.atoms:
                        self.atoms.add_attribute(
                            "formal_charge", coltype="int", default=0
                        )
                    ids = self.atoms.append(
                        x=Xs, y=Ys, z=Zs, atno=atnos, formal_charge=qs
                    )
                else:
                    ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)
            except Exception:
                print("Xs=")
                pprint.pp(Xs)
                print("Ys=")
                pprint.pp(Ys)
                print("Zs=")
                pprint.pp(Zs)
                print("atnos=")
                pprint.pp(atnos)
                raise
        else:
            ids = self.atoms.ids

            if coordinates:
                xyz = [[x, y, z] for x, y, z in zip(Xs, Ys, Zs)]
                self.atoms.coordinates = xyz

        if atoms or bonds:
            i = [ids[x] for x in Is]
            j = [ids[x] for x in Js]
            self.bonds.append(i=i, j=j, bondorder=BondOrders)

        data = mol.GetPropsAsDict()
        if self.__class__.__name__ == "_Configuration":
            # Charge
            self.charge = int(Chem.GetFormalCharge(mol))

            # Calculate spin multiplicity assuming maximal spin
            n_electrons = 0
            for rdk_atom in mol.GetAtoms():
                n_electrons += rdk_atom.GetNumRadicalElectrons()
            self.spin_multiplicity = n_electrons + 1

            # Check for property items for charge and multiplicity
            if "SEAMM|net charge|int|" in data:
                self.charge = int(data["SEAMM|net charge|int|"])
                del data["SEAMM|net charge|int|"]
            if "SEAMM|spin multiplicity|int|" in data:
                self.spin_multiplicity = int(data["SEAMM|spin multiplicity|int|"])
                del data["SEAMM|spin multiplicity|int|"]

            # Coordinates
            if coordinates and "SEAMM|XYZ|json|" in data:
                self.coordinates = json.loads(data["SEAMM|XYZ|json|"])
                del data["SEAMM|XYZ|json|"]

            if "SEAMM|system name|str|" in data:
                system_name = data["SEAMM|system name|str|"]
                del data["SEAMM|system name|str|"]
            if "SEAMM|configuration name|str|" in data:
                configuration_name = data["SEAMM|configuration name|str|"]
                del data["SEAMM|configuration name|str|"]

        # Record any properties in the database if desired
        if properties == "all":
            for key, value in data.items():
                if key.startswith("SEAMM|"):
                    # SEAMM|<property>|<type>|<units>. A property's name may hold
                    # '|' (density#LAMMPS#LAMMPS:MD|FF@oplsaa+), so split the
                    # type and units off the end.
                    _property, _type, units = key[len("SEAMM|") :].rsplit("|", 2)
                    # Keep dimensionless units as "", matching the property
                    # definitions, rather than turning them into None (NULL).
                    units = units.strip()
                    if not self.properties.exists(_property):
                        self.properties.add(_property, _type=_type, units=units)

                    if _type == "json":
                        value = json.dumps(value)

                    self.properties.put(_property, value)
                else:
                    if not self.properties.exists(key):
                        _type = value.__class__.__name__
                        self.properties.add(key, _type)

                    self.properties.put(key, value)
        return system_name, configuration_name

    def get_SSSR(self, symmetrized=False):
        """Calculate the smallest set of smallest rings.

        Parameters
        ----------
        symmetrized : bool = False
            Whether to return the symmetrized SSSR

        Returns
        -------
        list of list of int
        """
        if symmetrized:
            return Chem.rdmolops.GetSymSSSR(self.to_RDKMol())

        return Chem.rdmolops.GetSSSR(self.to_RDKMol())

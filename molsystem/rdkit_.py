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

# Valence
valence = {
    1: 1,
    2: 0,
    3: 1,
    4: 2,
    5: 3,
    6: 4,
    7: 3,
    8: 2,
    9: 1,
    10: 0,
    11: 1,
    12: 2,
    13: 3,
    14: 4,
    15: 3,
    16: 2,
    17: 1,
    18: 0,
    92: 2,
}

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
            rdk_mol.AddBond(
                index[row["i"]], index[row["j"]], bond_types[row["bondorder"]]
            )

        # Check for NH4+ type groups and set their charge
        rdk_mol.UpdatePropertyCache(strict=False)
        for at in rdk_mol.GetAtoms():
            atno = at.GetAtomicNum()
            if atno in valence:
                if atno == 5:
                    # Boron is a special case! Both BH4- and BH2- have a - charge.
                    if at.GetExplicitValence() == 4:
                        charge = -1
                    elif at.GetExplicitValence() == 2:
                        charge = -1
                    else:
                        charge = 0
                else:
                    charge = at.GetExplicitValence() - valence[atno]
                if charge != 0 and at.GetFormalCharge() == 0:
                    at.SetFormalCharge(charge)

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
        for rdk_atom in mol.GetAtoms():
            atnos.append(rdk_atom.GetAtomicNum())
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
            if "SEAMM|spin multiplicity" in data:
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
                    _, _property, _type, units = key.split("|", 4)
                    units = None if units.strip() == "" else units
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

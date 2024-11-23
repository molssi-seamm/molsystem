# -*- coding: utf-8 -*-

"""Interface to OpenEye OEChem."""

import logging

try:
    from openeye import oechem
except ImportError:
    oechem_available = False
    oechem_licensed = False
else:
    oechem_available = True
    oechem_licensed = oechem.OEChemIsLicensed()

logger = logging.getLogger(__name__)


def check_openeye_license():
    if not oechem_available:
        raise RuntimeError(
            "OpenEye OEChem is not installed! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )
    if not oechem_licensed:
        raise RuntimeError(
            "Cannot find a license for OpenEye OEChem! See "
            "https://docs.eyesopen.com/toolkits/python/index.html for detail"
        )


def openeye_version():
    """The version of the OpenEye OEChem toolkit."""
    if not oechem_available:
        return None
    return oechem.OEChemGetVersion()


class OpenEyeMixin:
    """A mixin for handling OpenEye's software via its Python interface."""

    def to_OEGraphMol(self, properties=None):
        """Return an OEGraphMol object for the configuration, template, or subset."""
        check_openeye_license()
        oe_mol = oechem.OEGraphMol()

        if self.__class__.__name__ == "_Configuration":
            oe_mol.SetIntData("net charge", self.charge)
            oe_mol.SetIntData("spin multiplicity", self.spin_multiplicity)

        # Create the atoms
        oe_atoms = []
        if "charge" in self.atoms:
            charges = self.atoms.get_column_data("charge")
            for atno, xyz, q in zip(
                self.atoms.atomic_numbers, self.atoms.coordinates, charges
            ):
                oe_atom = oe_mol.NewAtom(atno)
                oe_mol.SetCoords(oe_atom, xyz)
                oe_atom.SetPartialCharge(q)
                oe_atoms.append(oe_atom)
        else:
            for atno, xyz in zip(self.atoms.atomic_numbers, self.atoms.coordinates):
                oe_atom = oe_mol.NewAtom(atno)
                oe_mol.SetCoords(oe_atom, xyz)
                oe_atoms.append(oe_atom)

        # and bonds
        index = {j: i for i, j in enumerate(self.atoms.ids)}
        for row in self.bonds.bonds():
            oe_mol.NewBond(
                oe_atoms[index[row["i"]]], oe_atoms[index[row["j"]]], row["bondorder"]
            )

        if properties == "all":
            data = self.properties.get("all", include_system_properties=True)
            for key, value in data.items():
                _type = self.properties.type(key)
                if _type == "int":
                    oe_mol.SetIntData(key, value)
                elif _type == "float":
                    oe_mol.SetDoubleData(key, value)
                elif _type == "str":
                    oe_mol.SetStringData(key, value)
                else:
                    raise ValueError(f"Can't handle property of type '{_type}'")

                # Units, if any
                units = self.properties.units(key)
                if units is not None and units != "":
                    tmp = key.split("#", maxsplit=1)
                    if len(tmp) > 1:
                        oe_mol.SetStringData(tmp[0] + ",units" + "#" + tmp[1], units)
                    else:
                        oe_mol.SetStringData(key + ",units", units)
        return oe_mol

    def from_OEMol(
        self, oe_mol, properties="all", atoms=True, coordinates=True, bonds=True
    ):
        """Transform an OpenEye molecule into the current object."""
        check_openeye_license()

        atnos = []
        Xs = []
        Ys = []
        Zs = []
        qs = []
        for oe_atom in oe_mol.GetAtoms():
            atno = oe_atom.GetAtomicNum()
            atnos.append(atno)
            x, y, z = oe_mol.GetCoords(oe_atom)
            Xs.append(x)
            Ys.append(y)
            Zs.append(z)
            qs.append(oe_atom.GetPartialCharge())

        index = {at: i for i, at in enumerate(oe_mol.GetAtoms())}
        Is = []
        Js = []
        BondOrders = []
        for oe_bond in oe_mol.GetBonds():
            oe_i = oe_bond.GetBgn()
            oe_j = oe_bond.GetEnd()
            i = index[oe_i]
            j = index[oe_j]
            bondorder = oe_bond.GetOrder()
            Is.append(i)
            Js.append(j)
            BondOrders.append(bondorder)
            logger.debug(f"bond {i} - {j} {bondorder}")

        if atoms:
            self.clear()

        # Get the property data, cast to correct type
        data = {}
        for tmp in oe_mol.GetDataIter():
            tag = tmp.GetTag()
            attribute = oechem.OEGetTag(tag)
            value = oe_mol.GetData(tag)
            data[attribute] = value

        # Check for property items for charge and multiplicity
        if self.__class__.__name__ == "_Configuration":
            if "net charge" in data:
                self.charge = int(data["net charge"])
            if "spin multiplicity" in data:
                self.spin_multiplicity = int(data["spin multiplicity"])

        if atoms:
            if any([i != 0.0 for i in qs]):
                if "formal_charge" not in self.atoms:
                    self.atoms.add_attribute("charge", coltype="float", default=0)
                ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos, charge=qs)
            else:
                ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)
        else:
            ids = self.atoms.ids

            if coordinates:
                xyz = [[x, y, z] for x, y, z in zip(Xs, Ys, Zs)]
                self.atoms.coordinates = xyz

        if atoms or bonds:
            i = [ids[x] for x in Is]
            j = [ids[x] for x in Js]
            self.bonds.append(i=i, j=j, bondorder=BondOrders)

        # Record any properties in the database if desired
        if properties == "all":
            for key, value in data.items():
                if ",units" not in key and key not in [
                    "",
                    "net charge",
                    "spin multiplicity",
                ]:
                    if not self.properties.exists(key):
                        tmp = key.split("#", maxsplit=1)
                        if len(tmp) > 1:
                            units_key = tmp[0] + ",units" + "#" + tmp[1]
                        else:
                            units_key = key + ",units"
                        _type = value.__class__.__name__
                        if units_key in data:
                            units = data[units_key]
                            self.properties.add(key, _type, units=units)
                        else:
                            self.properties.add(key, _type)
                    self.properties.put(key, value)
        return self

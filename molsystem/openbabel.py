# -*- coding: utf-8 -*-

"""Interface to openbabel."""

import logging
from pathlib import Path

try:
    import openbabel  # noqa: F401
    import openbabel.openbabel as ob
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise

logger = logging.getLogger(__name__)


class OpenBabelMixin:
    """A mixin for handling OpenBabel via its Python interface."""

    def to_OBMol(self, properties=None):
        """Return an OBMol object for the configuration, template, or subset."""
        if self.__class__.__name__ == "_Configuration":
            charge = self.charge
            spin = self.spin_multiplicity
        else:
            charge = None
            spin = None

        ob_mol = ob.OBMol()
        for atno, xyz in zip(self.atoms.atomic_numbers, self.atoms.coordinates):
            ob_atom = ob_mol.NewAtom()
            ob_atom.SetAtomicNum(atno)
            ob_atom.SetVector(*xyz)
            if charge is not None:
                ob_atom.SetFormalCharge(charge)
                charge = None
            if spin is not None:
                ob_atom.SetSpinMultiplicity(spin)
                spin = None

        # 1-based indices in ob.
        index = {j: i for i, j in enumerate(self.atoms.ids, start=1)}
        for row in self.bonds.bonds():
            ob_mol.AddBond(index[row["i"]], index[row["j"]], row["bondorder"])

        if self.__class__.__name__ == "_Configuration":
            ob_mol.SetTotalCharge(self.charge)
            if self.spin_multiplicity is None:
                n_electrons = sum(self.atoms.atomic_numbers) - self.charge
                if n_electrons % 2 == 0:
                    multiplicity = 1
                else:
                    multiplicity = 2
                self.spin_multiplicity = multiplicity
            ob_mol.SetTotalSpinMultiplicity(self.spin_multiplicity)

        # Set local radical character
        ob_mol.AssignSpinMultiplicity(True)

        if properties == "all":
            data = self.properties.get("all", include_system_properties=True)
            pair = ob.OBPairData()
            for key, value in data.items():
                pair.SetAttribute(key)
                pair.SetValue(str(value))
                ob_mol.CloneData(pair)

        return ob_mol

    def from_OBMol(
        self, ob_mol, properties="all", atoms=True, coordinates=True, bonds=True
    ):
        """Transform an Open Babel molecule into the current object."""
        atnos = []
        Xs = []
        Ys = []
        Zs = []
        qs = []
        for ob_atom in ob.OBMolAtomIter(ob_mol):
            atno = ob_atom.GetAtomicNum()
            atnos.append(atno)
            Xs.append(ob_atom.x())
            Ys.append(ob_atom.y())
            Zs.append(ob_atom.z())
            logger.debug(f"atom {atno} {ob_atom.x()} {ob_atom.z()} {ob_atom.z()}")
            qs.append(ob_atom.GetFormalCharge())

        Is = []
        Js = []
        BondOrders = []
        for ob_bond in ob.OBMolBondIter(ob_mol):
            ob_i = ob_bond.GetBeginAtom()
            ob_j = ob_bond.GetEndAtom()
            i = ob_i.GetIdx()
            j = ob_j.GetIdx()
            bondorder = ob_bond.GetBondOrder()
            Is.append(i)
            Js.append(j)
            BondOrders.append(bondorder)
            logger.debug(f"bond {i} - {j} {bondorder}")

        if atoms:
            self.clear()

        if self.__class__.__name__ == "_Configuration":
            self.charge = ob_mol.GetTotalCharge()
            self.spin_multiplicity = ob_mol.GetTotalSpinMultiplicity()

        if atoms:
            if any([i != 0.0 for i in qs]):
                if "formal_charge" not in self.atoms:
                    self.atoms.add_attribute("formal_charge", coltype="int", default=0)
                ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos, formal_charge=qs)
            else:
                ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)
        else:
            ids = self.atoms.ids

            if coordinates:
                xyz = [[x, y, z] for x, y, z in zip(Xs, Ys, Zs)]
                self.atoms.coordinates = xyz

        if atoms or bonds:
            i = [ids[x - 1] for x in Is]
            j = [ids[x - 1] for x in Js]
            self.bonds.append(i=i, j=j, bondorder=BondOrders)

        # Record any properties in the database if desired
        if properties == "all":
            data = ob_mol.GetData()
            for item in data:
                attribute = item.GetAttribute()
                value = item.GetValue()
                if not self.properties.exists(attribute):
                    try:
                        int(value)
                        _type = "int"
                    except Exception:
                        try:
                            float(value)
                            _type = "float"
                        except Exception:
                            _type = "str"
                    self.properties.add(
                        attribute,
                        _type,
                        description="Imported from SDF file",
                    )
                self.properties.put(attribute, value)
        return self

    def coordinates_from_OBMol(self, ob_mol):
        """Update the coordinates from an Open Babel molecule."""
        XYZs = []
        for ob_atom in ob.OBMolAtomIter(ob_mol):
            XYZs.append([ob_atom.x(), ob_atom.y(), ob_atom.z()])

        self.coordinates = XYZs

    def coordinates_to_OBMol(self, ob_mol):
        """Update the coordinates of an Open Babel molecule from the configuration."""
        XYZs = self.coordinates
        for XYZ, ob_atom in zip(XYZs, ob.OBMolAtomIter(ob_mol)):
            ob_atom.SetVector(*XYZ)

    def find_substructures(self, template):
        """Find the substructures matching the template.

        Parameters
        ----------
        template : str, _Configuration, _Template, or _Subset
            The template, which may be a SMARTS string, or a molecular object.

        Returns
        -------
        [[int]]
            Lists of atom ids for matches.
        """
        if isinstance(template, str):
            smarts = template
        else:
            smarts = self.smiles

        ob_mol = self.to_OBMol()

        pattern = ob.OBSmartsPattern()
        pattern.Init(smarts)
        pattern.Match(ob_mol)
        maplist = pattern.GetUMapList()

        return [x for x in maplist]

    def from_sdf_text(self, text):
        """Get the text of an SDF file for the configuration.

        Parameters
        ----------
        text : str
            The text of an SDF file

        Returns
        -------
        (str, str)
            system name, configuration name
        """
        obConversion = ob.OBConversion()
        obConversion.SetInFormat("sdf")
        obMol = ob.OBMol()
        obConversion.ReadString(obMol, text)

        self.from_OBMol(obMol)

        # See if the system and configuration names are encoded in the tit
        result = (None, None)
        title = obMol.GetTitle()
        if "SEAMM=" in title:
            for tmp in title.split("|"):
                if "SEAMM=" in tmp and "/" in tmp:
                    sysname, confname = tmp.split("=", 1)[1].split("/")
                    sysname = sysname.strip()
                    confname = confname.strip()
                    result = (sysname, confname)
        return result

    def from_sdf(self, path):
        """Directly read an SDF file for the configuration.

        Parameters
        ----------
        path : pathlib.Path or str
            The path or name of the file to write.

        Returns
        -------
        (str, str)
            system name, configuration name
        """
        path = Path(path)
        text = path.read_text()

        return self.from_sdf_text(text)

    def to_sdf_text(self):
        """Get the text of an SDF file for the configuration.

        Returns
        -------
        str
            The text of the SDF file
        """
        obConversion = ob.OBConversion()
        obConversion.SetOutFormat("sdf")
        obMol = self.to_OBMol(properties="all")
        title = f"SEAMM={self.system.name}/{self.name}"
        obMol.SetTitle(title)
        text = obConversion.WriteString(obMol)

        return text

    def to_sdf(self, path):
        """Directly write an SDF file for the configuration.

        Parameters
        ----------
        path : pathlib.Path or str
            The path or name of the file to write.

        Returns
        -------
        str
            The text of the SDF file
        """
        text = self.to_sdf_text()

        path = Path(path)
        path.write_text(text)

        return text

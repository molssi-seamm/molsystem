# -*- coding: utf-8 -*-

"""Interface to openbabel."""

import logging

try:
    from openbabel import openbabel
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
        ob_mol = openbabel.OBMol()
        for atno, xyz in zip(self.atoms.atomic_numbers, self.atoms.coordinates):
            ob_atom = ob_mol.NewAtom()
            ob_atom.SetAtomicNum(atno)
            ob_atom.SetVector(*xyz)

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
            pair = openbabel.OBPairData()
            for key, value in data.items():
                pair.SetAttribute(key)
                pair.SetValue(str(value))
                ob_mol.CloneData(pair)

        return ob_mol

    def from_OBMol(self, ob_mol, properties="all"):
        """Transform an Open Babel molecule into the current object."""
        atnos = []
        Xs = []
        Ys = []
        Zs = []
        for ob_atom in openbabel.OBMolAtomIter(ob_mol):
            atno = ob_atom.GetAtomicNum()
            atnos.append(atno)
            Xs.append(ob_atom.x())
            Ys.append(ob_atom.y())
            Zs.append(ob_atom.z())
            logger.debug(f"atom {atno} {ob_atom.x()} {ob_atom.z()} {ob_atom.z()}")

        Is = []
        Js = []
        BondOrders = []
        for ob_bond in openbabel.OBMolBondIter(ob_mol):
            ob_i = ob_bond.GetBeginAtom()
            ob_j = ob_bond.GetEndAtom()
            i = ob_i.GetIdx()
            j = ob_j.GetIdx()
            bondorder = ob_bond.GetBondOrder()
            Is.append(i)
            Js.append(j)
            BondOrders.append(bondorder)
            logger.debug(f"bond {i} - {j} {bondorder}")

        self.clear()

        if self.__class__.__name__ == "_Configuration":
            self.charge = ob_mol.GetTotalCharge()
            self.spin_multiplicity = ob_mol.GetTotalSpinMultiplicity()

        ids = self.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)
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

        pattern = openbabel.OBSmartsPattern()
        pattern.Init(smarts)
        pattern.Match(ob_mol)
        maplist = pattern.GetUMapList()

        return [x for x in maplist]

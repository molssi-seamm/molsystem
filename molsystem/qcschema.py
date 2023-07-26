# -*- coding: utf-8 -*-

"""Interface to qcschema."""

import json
import logging

from seamm_util import Q_

logger = logging.getLogger(__name__)


class QCSchemaMixin:
    """A mixin for handling QCSchema."""

    def to_qcschema_dict(self, properties=None):
        """Create a dictionary compliant with QCSchema."""
        result = {
            "schema_name": "qcschema_molecule",
            "schema_version": 2,
        }

        # Symbols and coordinates (in Bohr)
        result["symbols"] = [*self.atoms.symbols]
        factor = Q_(1.0, "Å").m_as("a_0")
        xyz = []
        # round() below helps tests work across platforms. 9 digits are enough!
        for row in self.atoms.get_coordinates(fractionals=False):
            for val in row:
                xyz.append(round(val * factor, 9))
        result["geometry"] = xyz

        # Charge and multiplicity
        result["molecular_charge"] = self.charge
        result["molecular_multiplicity"] = self.spin_multiplicity

        # Bonds, if any
        bonds = []
        index = {j: i for i, j in zip(range(self.n_atoms), self.atoms.ids)}
        for row in self.bonds.bonds():
            bonds.append((index[row["i"]], index[row["j"]], row["bondorder"]))
        if len(bonds) > 0:
            result["connectivity"] = bonds

        # Molecules (fragments in QCSchema speak)
        result["fragments"] = self.find_molecules(as_indices=True)

        result["name"] = f"{self.system.name} / {self.name}"
        if "name" in self.atoms:
            result["atom_labels"] = self.atoms.get_column_data("name")

        return result

    def to_qcschema_json(self):
        """Create the QCSchema JSON for the molecule."""
        data = self.to_qcschema_dict()
        return json.dumps(data)

    def from_qcschema_dict(self, data):
        """Reset the molecule from the QCSchema data."""
        self.clear()
        self.periodicity = 0

        symbols = data["symbols"]
        factor = Q_(1.0, "a_0").m_as("Å")
        Xs = []
        Ys = []
        Zs = []
        for x, y, z in zip(*[iter(data["geometry"])] * 3):
            Xs.append(x * factor)
            Ys.append(y * factor)
            Zs.append(z * factor)

        ids = self.atoms.append(x=Xs, y=Ys, z=Zs, symbol=symbols)

        if "atom_labels" in data:
            if "name" not in self.atoms:
                self.atoms.add_attribute("name", values=data["atom_labels"])
            else:
                self.atoms["name"] = data["atom_labels"]

        if "connectivity" in data and len(data["connectivity"]) > 0:
            Is = []
            Js = []
            orders = []
            for i, j, order in data["connectivity"]:
                Is.append(i)
                Js.append(j)
                orders.append(order)

            i = [ids[x - 1] for x in Is]
            j = [ids[x - 1] for x in Js]
            self.bonds.append(i=i, j=j, bondorder=orders)

    def from_qcschema_json(self, json_data):
        """Reset the molecule from the QCSchema JSON."""
        data = json.loads(json_data)
        self.from_qcschema_dict(data)

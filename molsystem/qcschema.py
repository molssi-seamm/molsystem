# -*- coding: utf-8 -*-

"""Interface to qcschema."""

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
        for row in self.atoms.get_coordinates(fractionals=False):
            for val in row:
                xyz.append(val * factor)
        result["geometry"] = xyz

        # Charge and multiplicity
        result["molecular_charge"] = self.charge
        result["molecular_multiplicity"] = self.spin_multiplicity

        # Bonds, if any
        bonds = []
        for row in self.bonds.bonds():
            bonds.append((row["i"], row["j"], row["bondorder"]))
        if len(bonds) > 0:
            result["connectivity"] = bonds

        # Molecules (fragments in QCSchema speak)
        result["fragments"] = self.find_molecules(as_indices=True)

        result["name"] = self.name
        if "name" in self.atoms:
            result["atom_labels"] = self.atoms.get_column_data("name")

        return result

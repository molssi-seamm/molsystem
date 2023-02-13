# -*- coding: utf-8 -*-

"""Interface to the CMS schema."""

import logging

logger = logging.getLogger(__name__)


class CMSSchemaMixin:
    """A mixin for handling CMS Schema."""

    def to_cms_schema(self, properties=None):
        """Create a dictionary compliant with CMS Schema."""
        data = {}

        # Symbols and coordinates
        data["symbols"] = [*self.atoms.symbols]
        data["periodicity"] = self.periodicity
        data["charge"] = self.charge
        data["multiplicity"] = self.spin_multiplicity
        if self.periodicity == 3:
            data["cell"] = self.cell.parameters
            data["spacegroup"] = self.symmetry.group

        coordinates = data["coordinates"] = {}
        coordinates["coordinate system"] = self.coordinate_system
        coordinates["units"] = "Å"
        coordinates["coordinates"] = self.atoms.get_coordinates()

        # Bonds, if any
        bonds = []
        for row in self.bonds.bonds():
            bonds.append((row["i"], row["j"], row["bondorder"]))
        if len(bonds) > 0:
            data["bonds"] = bonds

        data["name"] = self.name
        if "name" in self.atoms:
            data["atom labels"] = self.atoms.get_column_data("name")

        result = {
            "schema name": "cms_schema",
            "schema version": "1.0",
            "systems": [
                {
                    "name": self.system.name,
                    "configurations": [data],
                },
            ],
        }

        return result

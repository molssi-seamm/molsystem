# -*- coding: utf-8 -*-

"""Functions for handling InChI"""

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


class InChIMixin:
    """A mixin for handling InChI."""

    @property
    def inchi(self):
        """Return the InChI string for this object."""
        return self.to_inchi()

    @property
    def inchikey(self):
        """Return the InChIKey string for this object."""
        return self.to_inchi(key=True)

    def to_inchi(self, key=False):
        """Create the InChI string from the system.

        Parameters
        ----------
        key : bool = False
            Whether to create the InChIKey

        Returns
        -------
        str
            The InChI string, or (InChI, name) if the rname is requested
        """
        logger.info("to_inchi")

        obConversion = openbabel.OBConversion()
        if key:
            obConversion.SetOutFormat("inchikey")
        else:
            obConversion.SetOutFormat("inchi")

        mol = self.to_OBMol()

        inchi = obConversion.WriteString(mol)

        logger.info(f"inchi = '{inchi}'")

        return inchi.strip()

    def from_inchi(self, inchi, name=None):
        """Create the system from a InChI string.

        Parameters
        ----------
        inchi : str
            The InChI string
        name : str = None
            The name of the molecule

        Returns
        -------
        None
        """

        save = self.name

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("inchi", "mdl")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, inchi)

        # Add hydrogens
        mol.AddHydrogens()

        # Get coordinates for a 3-D structure
        builder = openbabel.OBBuilder()
        builder.Build(mol)

        self.from_OBMol(mol)

        if name is not None:
            self.name = name
        else:
            self.name = save

# -*- coding: utf-8 -*-

"""Functions for handling InChI"""

import logging
import requests

try:
    from openbabel import openbabel
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "     conda install -c conda-forge openbabel"
    )
    raise
try:
    from rdkit import Chem
except ModuleNotFoundError:
    print("Please install RDKit using conda:\n     conda install -c conda-forge rdkit")
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

    def to_inchi(self, key=False, rdkit=False):
        """Create the InChI string from the system.

        Parameters
        ----------
        key : bool = False
            Whether to create the InChIKey
        rdkit : bool = False
            Whether to use RDKit rather than default of OpenBabel

        Returns
        -------
        str
            The InChI string, or (InChI, name) if the rname is requested
        """
        logger.info("to_inchi")

        if rdkit:
            mol = self.to_RDKMol()
            if key:
                inchi = Chem.inchi.MolToInchiKey(mol)
            else:
                inchi = Chem.inchi.MolToInchi(mol)
        else:
            obConversion = openbabel.OBConversion()
            if key:
                obConversion.SetOutFormat("inchikey")
            else:
                obConversion.SetOutFormat("inchi")

            mol = self.to_OBMol()

            inchi = obConversion.WriteString(mol)

        logger.info(f"inchi = '{inchi}'")

        return inchi.strip()

    def from_inchi(self, inchi, name=None, rdkit=False):
        """Create the system from a InChI string.

        Parameters
        ----------
        inchi : str
            The InChI string
        name : str = None
            The name of the molecule
        rdkit : bool = False
            Whether to use RDKit rather than default of OpenBabel

        Returns
        -------
        None
        """

        save = self.name

        if rdkit:
            mol = Chem.inchi.MolFromInchi(inchi, sanitize=True, removeHs=False)
            self.from_RDKMol(mol)
        else:
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

    def from_inchikey(self, inchikey, name=None):
        """Create the system from an InChIKey string.

        Parameters
        ----------
        inchikey : str
            The InChIKey string
        name : str = None
            The name of the molecule

        Returns
        -------
        None
        """
        inchi = self._get_inchi(inchikey)

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

    def _get_inchi(self, inchikey):
        """Get the InChI from PubChem given the InChIKey."""
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}"
            "/property/InChI/JSON"
        )
        r = requests.get(url)
        result = r.json()
        if "Fault" in result:
            raise RuntimeError(f"InChIKey '{inchikey}' not found in PubChem.")

        inchi = set()
        for properties in result["PropertyTable"]["Properties"]:
            if "InChI" in properties:
                inchi.add(properties["InChI"])
        if len(inchi) == 1:
            return inchi.pop()
        else:
            return [*inchi]

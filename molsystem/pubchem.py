# -*- coding: utf-8 -*-

"""Functions for handling PubChem"""

import logging
from urllib.parse import quote as url_quote

import requests

logger = logging.getLogger(__name__)
pug_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

try:
    import pubchempy as pcp
except ModuleNotFoundError:
    print(
        "To use PubChem, please install pubchempy using conda:\n"
        "     conda install -c conda-forge pubchempy"
    )
    logger.warning(
        "To use PubChem, please install pubchempy using conda:\n"
        "     conda install -c conda-forge pubchempy"
    )


class PubChemMixin:
    """A mixin for handling the PubChem database."""

    @property
    def PC_cid(self):
        """Return the PubChem CID for this structure, or None."""
        inchi = self.to_inchi()
        results = pcp.get_compounds(inchi, "inchi")

        if len(results) == 0:
            return None

        if len(results) >= 1:
            logger.info(f"PubChem search returned more than one hit of {inchi}")

        result = results[0].to_dict()
        if "cid" in result:
            return result["cid"]

        return None

    def PC_from_cid(self, cid, fallback=None):
        """Create the configuration from the PubChem 3-D structure, if available.

        Parameters
        ----------
        cid : int
            The PubChem CID.

        fallback : str
            A fallback SMILES, InChI, etc. to use if PubChem fails
        """
        response = requests.get(f"{pug_url}/compound/cid/{cid}/SDF")
        if response.status_code == 200:
            self.from_sdf_text(response.text)
            return

        # An error!
        if fallback is None:
            raise RuntimeError(f"No 3-D structure available for {cid}")

        # See what the fallback is
        if len(fallback) == 27:
            tmp = fallback.split("-")
            if (
                len(tmp) == 3
                and len(tmp[0]) == 14
                and len(tmp[1]) == 10
                and len(tmp[2]) == 1
            ):
                self.from_inchikey(fallback)
                return
        if fallback[0:7] == "InChI=":
            self.from_inchi(fallback)
            return
        else:
            self.from_smiles(fallback)
            return

    def PC_from_identifier(self, identifier, namespace="detect", fallback=None):
        """Create the configuration from the PubChem 3-D structure, if available.

        Parameters
        ----------
        identifier : int or str
            The PubChem identifier
        namespace : str
            The PubChem namespace: cid, name, smiles, inchi, inchikey
        fallback : str
            A fallback SMILES, InChI, etc. to use if PubChem fails
        """
        if namespace == "detect":
            # Work through the possibilities
            if len(identifier) == 27:
                tmp = identifier.split("-")
                if (
                    len(tmp) == 3
                    and len(tmp[0]) == 14
                    and len(tmp[1]) == 10
                    and len(tmp[2]) == 1
                ):
                    namespace = ["inchikey"]
            elif identifier[0:7] == "InChI=":
                self.from_inchi(identifier)
                namespaces = ["inchi"]
            else:
                namespaces = ["name", "smiles"]
        else:
            namespaces = [namespace]

        for namespace in namespaces:
            response = requests.get(
                f"{pug_url}/compound/{namespace}/{url_quote(identifier)}/SDF"
            )
            if response.status_code == 200:
                self.from_sdf_text(response.text)
                return

        # An error!
        if fallback is None:
            raise RuntimeError(f"No 3-D structure available for {identifier}")

        # See what the fallback is
        if len(fallback) == 27:
            tmp = fallback.split("-")
            if (
                len(tmp) == 3
                and len(tmp[0]) == 14
                and len(tmp[1]) == 10
                and len(tmp[2]) == 1
            ):
                self.from_inchikey(fallback)
                return
        if fallback[0:7] == "InChI=":
            self.from_inchi(fallback)
            return
        else:
            self.from_smiles(fallback)
            return

    def PC_iupac_name(self, fallback=None):
        """Return the IUPAC name for this structure, or None.

        Parameters
        ----------
        fallback : str
            A name to return if PubChem doesn't have a name
        """
        inchi = self.to_inchi()
        results = pcp.get_compounds(inchi, "inchi")

        if len(results) == 0:
            return fallback

        if len(results) >= 1:
            logger.info(f"PubChem search returned more than one hit of {inchi}")

        result = results[0].to_dict()
        if "iupac_name" in result:
            return result["iupac_name"]

        return fallback

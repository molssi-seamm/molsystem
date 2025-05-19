# -*- coding: utf-8 -*-

"""Functions for handling PubChem"""

import logging
import re
from time import perf_counter_ns, sleep
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


def PC_standardize(structures):
    """Use the PubChem standardization service with SMILES.

    Parameters
    ----------
    structures : [str]
        The SMILES of the structures to standardize

    Returns
    -------
    [{str: str}]
        A list of dicts of the results.

    Note:
        The response contains information to throttle the requests, like this:

        X-Throttling-Control: Request Count status: Green (0%),
                              Request Time status: Green (0%),
                              Service status: Green (20%)
    """
    results = []
    max_count = 0
    max_time = 0
    max_service = 0
    n_sleep = 0
    with requests.Session() as session:
        t0 = perf_counter_ns()
        for SMILES in structures:
            while True:
                response = session.post(
                    f"{pug_url}/standardize/SMILES/JSON",
                    data={"smiles": SMILES, "include_components": False},
                )

                # Check the throttling request
                if "X-Throttling-Control" in response.headers:
                    throttling_header = response.headers["X-Throttling-Control"]
                    match = re.search(
                        r"Request Count status: \w+ \((\d+)%\), "
                        r"Request Time status: \w+ \((\d+)%\), "
                        r"Service status: \w+ \((\d+)%\)",
                        throttling_header,
                    )
                    if match:
                        count = int(match.group(1))
                        time = int(match.group(2))
                        service = int(match.group(3))

                        max_count = max(count, max_count)
                        max_time = max(time, max_time)
                        max_service = max(service, max_service)

                        if count > 75 or time > 75 or service > 75:
                            n_sleep += 1
                            sleep(1)

                if response.status_code == 503:
                    print("status code = 503")
                    continue

                response.raise_for_status()

                break

            data = response.json()

            # Pull out some info from the data
            cmpds = data["PC_Compounds"]
            if len(cmpds) > 1:
                raise ValueError(
                    f"There are {len(cmpds)} in the PubChem standardization!"
                )

            cmpd = cmpds[0]
            result = {"original SMILES": SMILES}
            if "id" in cmpd and "id" in cmpd["id"] and "cid" in cmpd["id"]["id"]:
                result["cid"] = cmpd["id"]["id"]["cid"]

            for props in cmpd["props"]:
                if "urn" in props and "value" in props and "label" in props["urn"]:
                    key = props["urn"]["label"]
                    value = props["value"]
                    _type = [*value.keys()][0]
                    if _type == "sval":
                        result[key] = value[_type]
                    else:
                        raise TypeError(f"Unknown type {_type} in compound properties.")

            results.append(result)
        t1 = perf_counter_ns()
        t = (t1 - t0) / 1000000000
        n = len(structures)
        per = t / n

    return {
        "data": results,
        "n_structures": len(structures),
        "throttling": {
            "maximum count": max_count,
            "maximum time": max_time,
            "maximum service": max_service,
            "number of slowdowns": n_sleep,
        },
        "time": {
            "total": t,
            "per structure": per,
        },
    }


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

    def PC_from_identifier(
        self, identifier, namespace="detect", fallback=None, properties="all"
    ):
        """Create the configuration from the PubChem 3-D structure, if available.

        Parameters
        ----------
        identifier : int or str
            The PubChem identifier
        namespace : str
            The PubChem namespace: cid, name, smiles, inchi, inchikey
        fallback : str
            A fallback SMILES, InChI, etc. to use if PubChem fails
        properties : str = "all"
            Whether to include all properties or none
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
                f"{pug_url}/compound/{namespace}/{url_quote(identifier)}/"
                "SDF?record_type=3d"
            )
            if response.status_code == 200:
                self.from_sdf_text(response.text, properties=properties)
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

# -*- coding: utf-8 -*-

"""Alignment methods for configurations"""

import logging
from math import sqrt

import numpy as np

try:
    import openbabel  # noqa: F401
    import openbabel.openbabel as ob
except ModuleNotFoundError:
    print(
        "Please install openbabel using conda:\n"
        "conda install -c conda-forge openbabel"
    )
    raise

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign, rdmolops, AllChem
except ModuleNotFoundError:
    print("Please install rdkit using conda:\n     conda install -c conda-forge rdkit")
    raise

logger = logging.getLogger(__name__)


def RMSD(
    structure, reference, include_h=False, symmetry=False, align=False, flavor=None
):
    """Compute the RMSD between configurations.

    Parameters
    ----------
    structure : Configuration, RDKMol, or OBMol

    reference : Configuration, RDKMol, or OBMol
        A single configuration to match. if 'structure' is not a configuration,
        this must match the type of 'structure'.

    include_h : bool = True
        Whether to include hydrogen atoms in the RMSD

    symmetry : bool = False
        Whether to detect symmetric flips. Note this requires a lot of memory for
        larger systems.

    align : bool = False
        Whether to transform the structure to the best alignment with
        the reference.

    flavor : str = None
        The flavor to use. Currently 'rdkit' or 'openbabel'

    Returns
    -------
    {str: float|int}
        RMSD: The RMS between the structure and the reference.
        maximum displacement: The maximum displacement of any atom.
        displaced atom: The index of the atom with the maximum displacement.
    """
    from .configuration import _Configuration

    if isinstance(structure, _Configuration):
        if isinstance(reference, _Configuration):
            if flavor is None or flavor == "rdkit":
                _structure = structure.to_RDKMol()
                _reference = reference.to_RDKMol()
            elif flavor == "openbabel":
                _structure = structure.to_OBMol()
                _reference = reference.to_OBMol()
            else:
                raise ValueError(f"Cannot handle flavor '{flavor}'")
        elif isinstance(reference, Chem.RWMol):
            if flavor is None or flavor == "rdkit":
                _structure = structure.to_RDKMol()
                _reference = Chem.RWMol(reference)
            else:
                raise ValueError(
                    f"Cannot use flavor '{flavor}' with an RDKMol reference"
                )
        elif isinstance(reference, ob.OBMol):
            if flavor is None or flavor == "openbabel":
                _structure = structure.to_OBMol()
                _reference = ob.OBMol(reference)
            else:
                raise ValueError(
                    f"Cannot use flavor '{flavor}' with an OBMol reference"
                )
        else:
            raise ValueError(
                "Reference must be a Configuration, RDKMol, or OBMol, not "
                f"{type(reference)}"
            )
    elif isinstance(structure, rdkit.Chem.rdchem.RWMol):
        if flavor is not None and flavor != "rdkit":
            raise ValueError(f"Cannot use flavor '{flavor}' with an RDKMol structure")
        if isinstance(reference, _Configuration):
            _structure = Chem.rdchem.RWMol(structure)
            _reference = reference.to_RDKMol()
        elif isinstance(reference, rdkit.Chem.rdchem.RWMol):
            _structure = Chem.rdchem.RWMol(structure)
            _reference = Chem.rdchem.RWMol(reference)
        elif isinstance(reference, ob.OBMol):
            raise ValueError(
                "If the structure is an RDKMol, the reference "
                f"must be a Configuration or RDKMol, not {type(reference)}"
            )
        else:
            raise ValueError(
                "Reference must be a Configuration, RDKMol, or OBMol, not "
                f"{type(reference)}"
            )
    elif isinstance(structure, ob.OBMol):
        if flavor is not None and flavor != "openbabel":
            raise ValueError(f"Cannot use flavor '{flavor}' with an OBMol structure")
        if isinstance(reference, _Configuration):
            _structure = ob.OBMol(structure)
            _reference = reference.to_OBMol()
        elif isinstance(reference, rdkit.Chem.rdchem.RWMol):
            raise ValueError(
                "If the structure is an OBMol, the reference "
                f"must be a Configuration or OBMol, not {type(reference)}"
            )
        elif isinstance(reference, ob.OBMol):
            _structure = ob.OBMol(structure)
            _reference = ob.OBMol(reference)
        else:
            raise ValueError(
                f"Reference must be a Configuration, or OBMol, not {type(reference)}"
            )
    else:
        raise ValueError(
            "Structure must be a Configuration, RDKMol, or OBMol, "
            f"not {type(structure)}"
        )

    if isinstance(_structure, rdkit.Chem.rdchem.RWMol):
        result = RDK_RMSD(
            _structure,
            _reference,
            include_h=include_h,
            symmetry=symmetry,
        )
        if align:
            if isinstance(structure, _Configuration):
                structure.coordinates_from_RDKMol(_structure)
            else:
                structure.GetConformer(0).SetPositions(
                    _structure.GetConformer(0).GetPositions()
                )
    elif isinstance(_structure, ob.OBMol):
        result = OB_RMSD(
            _structure,
            _reference,
            include_h=include_h,
            symmetry=symmetry,
        )
        if align:
            if isinstance(structure, _Configuration):
                structure.coordinates_from_OBMol(_structure)
            else:
                for at1, at2 in zip(
                    ob.OBMolAtomIter(structure), ob.OBMolAtomIter(_structure)
                ):
                    at1.SetVector(at2.GetVector())

    return result


def RDK_RMSD(structure, reference, include_h=True, symmetry=False):
    """Compute the RMSD between configurations using RDKit.

    Parameters
    ----------
    structure: RDKMol
        The structure to align

    reference : RDKMol
        The structure to compare to.

    include_h : bool = True
        Whether to include hydrogen atoms in the RMSD

    symmetry : bool = False
        Whether to detect symmetric flips. Note this requires a lot of memory for
        larger systems.

    Returns
    -------
    {str: float|int}
        RMSD: The RMS between the structure and the reference.
        maximum displacement: The maximum displacement of any atom.
        displaced atom: The index of the atom with the maximum displacement.

    Note
    ----
    The structure is modified in place to be optimally aligned with the
    reference.
    """
    if not include_h:
        structure = rdmolops.RemoveHs(structure)
        reference = rdmolops.RemoveHs(reference)

    if symmetry:
        rmsd, transform, atom_map = rdMolAlign.GetBestAlignmentTransform(
            structure, reference
        )
        AllChem.TransformMol(structure, transform)

        max_displacement = -1.0
        displaced_atom = None
        XYZ1 = structure.GetConformer(0).GetPositions()
        XYZ2 = reference.GetConformer(0).GetPositions()
        for at1, at2 in atom_map:
            dist = np.linalg.norm(XYZ1[at1] - XYZ2[at2])
            if dist > max_displacement:
                max_displacement = float(dist)
                displaced_atom = at1
    else:
        rmsd, transform = rdMolAlign.GetAlignmentTransform(structure, reference)
        AllChem.TransformMol(structure, transform)

        max_displacement = 0.0
        displaced_atom = None
        XYZ1 = structure.GetConformer(0).GetPositions()
        XYZ2 = reference.GetConformer(0).GetPositions()
        for count, XYZ1_XYZ2 in enumerate(zip(XYZ1, XYZ2)):
            XYZ1, XYZ2 = XYZ1_XYZ2
            dist = np.linalg.norm(XYZ1 - XYZ2)
            if dist > max_displacement:
                max_displacement = float(dist)
                displaced_atom = count
    return {
        "RMSD": rmsd,
        "maximum displacement": max_displacement,
        "displaced atom": displaced_atom,
    }


def OB_RMSD(structure, reference, include_h=True, symmetry=False):
    """Compute the RMSD between configurations using OpenBabel.

    Parameters
    ----------
    structure: OBMol
        The structure to align

    reference : OBMol
        The structure to compare to.

    include_h : bool = True
        Whether to include hydrogen atoms in the RMSD

    symmetry : bool = False
        Whether to detect symmetric flips. Note this requires a lot of memory for
        larger systems.

    Returns
    -------
    {str: float|int}
        RMSD: The RMS between the structure and the reference.
        maximum displacement: The maximum displacement of any atom.
        displaced atom: The index of the atom with the maximum displacement.

    Note
    ----
    The structure is modified in place to be optimally aligned with the
    reference.
    """
    align = ob.OBAlign(include_h, symmetry)
    align.SetTargetMol(structure)

    try:
        align.SetRefMol(reference)
        align.Align()
        align.UpdateCoords(structure)
        rmsd = align.GetRMSD()

        max_displacement = -1.0
        displaced_atom = None
        for count, at1_at2 in enumerate(
            zip(ob.OBMolAtomIter(structure), ob.OBMolAtomIter(reference))
        ):
            at1, at2 = at1_at2
            dist = sqrt(
                (at1.x() - at2.x()) ** 2
                + (at1.y() - at2.y()) ** 2
                + (at1.z() - at2.z()) ** 2
            )
            if dist > max_displacement:
                max_displacement = float(dist)
                displaced_atom = count
    except TypeError:
        raise
    except Exception as e:
        logger.error(f"Configuration.RMSD {e}")
        raise
    return {
        "RMSD": rmsd,
        "maximum displacement": max_displacement,
        "displaced atom": displaced_atom,
    }


class AlignMixin:
    """A mixin for handling alignment of configuration."""

    def RMSD(
        self, reference, include_h=False, symmetry=False, align=False, flavor="rdkit"
    ):
        """Compute the RMSD between configurations.

        Parameters
        ----------
        reference : Configuration|OBMol|RDKMol
            The structure to compare and possibly align to

        include_h : bool = True
            Whether to include hydrogen atoms in the RMSD

        symmetry : bool = False
            Whether to detect symmetric flips. Note this requires a lot of memory for
            larger systems.

        align : bool = False
            Whether to transform the structure to the best alignment with
            the target.

        flavor : str = "rdkit"
            The flavor to use. Currently 'rdkit' or 'openbabel'

        Returns
        -------
        float or [float]
            The RMS between the structure configuration and the target(s).
        """
        return RMSD(
            self,
            reference,
            include_h=include_h,
            symmetry=symmetry,
            align=align,
            flavor=flavor,
        )

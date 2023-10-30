# -*- coding: utf-8 -*-

"""Alignment methods for configurations"""

import logging

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


class AlignMixin:
    """A mixin for handling alignment of configuration."""

    def RMSD(self, other, include_h=True, symmetry=False):
        """Compute the RMSD between configurations.

        Parameters
        ----------
        other : Configuration or iterable of Configurations
            A single configuration or e.g. list of configurations to match

        include_h : bool = True
            Whether to include hydrogen atoms in the RMSD

        symmetry : bool = False
            Whether to detect symmetric flips. Note this requires a lot of memory for
            larger systems.

        Returns
        -------
        float or [float]
            The RMS between the current configuration and the target(s).
        """
        from .configuration import _Configuration

        align = ob.OBAlign(include_h, symmetry)
        align.SetRefMol(self.to_OBMol())

        if isinstance(other, _Configuration):
            align.SetTargetMol(other.to_OBMol())
            align.Align()
            return align.GetRMSD()
        elif isinstance(other, ob.OBMol):
            align.SetTargetMol(other)
            align.Align()
            return align.GetRMSD()

        result = []
        try:
            for target in other:
                if isinstance(target, _Configuration):
                    align.SetTargetMol(target.to_OBMol())
                    align.Align()
                    result.append(align.GetRMSD())
                elif isinstance(target, ob.OBMol):
                    align.SetTargetMol(target)
                    align.Align()
                    result.append(align.GetRMSD())
                else:
                    raise TypeError(f"RMSD can't handle {type(target)}.")
        except TypeError:
            raise
        except Exception as e:
            logger.error(f"Configuration.RMSD {e}")
            raise
        return result

    def align(self, other, include_h=True, symmetry=False):
        """Align configurations.

        Parameters
        ----------
        other : Configuration or iterable of Configurations
            A single configuration or e.g. list of configurations to match

        include_h : bool = True
            Whether to include hydrogen atoms in the RMSD

        symmetry : bool = False
            Whether to detect symmetric flips. Note this requires a lot of memory for
            larger systems.

        Returns
        -------
        float or [float]
            The RMS between the current configuration and the target(s).
        """
        from .configuration import _Configuration

        align = ob.OBAlign(include_h, symmetry)
        align.SetRefMol(self.to_OBMol())

        if isinstance(other, _Configuration):
            mol = other.to_OBMol()
            align.SetTargetMol(mol)
            align.Align()
            align.UpdateCoords(mol)
            other.coordinates_from_OBMol(mol)
            return align.GetRMSD()
        elif isinstance(other, ob.OBMol):
            align.SetTargetMol(other)
            align.Align()
            align.UpdateCoords(other)
            return align.GetRMSD()

        result = []
        try:
            for target in other:
                if isinstance(target, _Configuration):
                    mol = target.to_OBMol()
                    align.SetTargetMol(mol)
                    align.Align()
                    align.UpdateCoords(mol)
                    target.coordinates_from_OBMol(mol)
                    result.append(align.GetRMSD())
                elif isinstance(target, ob.OBMol):
                    align.SetTargetMol(target)
                    align.Align()
                    align.UpdateCoords(target)
                    result.append(align.GetRMSD())
                else:
                    raise TypeError(f"RMSD can't handle {type(target)}.")
        except TypeError:
            raise
        except Exception as e:
            logger.error(f"Configuration.RMSD {e}")
            raise
        return result

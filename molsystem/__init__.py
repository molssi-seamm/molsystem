# -*- coding: utf-8 -*-

"""
molsystem
A general implementation for molecular and periodic systems.
"""

# Temporary fix for issue with openeye and openbabel conflicting
# Seems to work if we import openeye first
# May not have openeye or openbabel installed, so take care
try:
    from openeye import oechem  # noqa: F401
except ImportError:
    pass
try:
    from openbabel import openbabel as ob  # noqa: F401
except ImportError:
    pass

# Bring up the classes so that they appear to be directly in
# the molsystem package.

import molsystem.elements  # noqa: F401
from .system_db import SystemDB  # noqa: F401
from .configuration import spin_states  # noqa: F401
from .cell import Cell  # noqa: F401
from .properties import standard_properties, add_properties_from_file  # noqa: F401
from .openbabel import openbabel_version, openbabel_citations  # noqa: F401
from .align import RMSD  # noqa: F401
from .pubchem import PC_standardize  # noqa: F401

try:
    from .openeye import openeye_version  # noqa: F401
except ImportError:
    pass
from .rdkit_ import rdkit_version, rdkit_citations  # noqa: F401


# Handle versioneer
from ._version import get_versions

__author__ = """Paul Saxe"""
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

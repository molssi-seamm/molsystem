# -*- coding: utf-8 -*-

"""Top-level package for MolSystem."""

__author__ = """Paul Saxe"""
__email__ = 'psaxe@vt.edu'
__version__ = '0.1.2'

# Bring up the classes so that they appear to be directly in
# the molsystem package.

from molsystem.table import Table  # noqa: F401
from molsystem.atoms import Atoms  # noqa: F401
from molsystem.bonds import Bonds  # noqa: F401
from molsystem.cell import Cell  # noqa: F401
from molsystem.system import System  # noqa: F401

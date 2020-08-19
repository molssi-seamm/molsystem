# -*- coding: utf-8 -*-

"""
molsystem
A general implementation for molecular and periodic systems.
"""

# Bring up the classes so that they appear to be directly in
# the molsystem package.

from molsystem.system import System  # noqa: F401

# Handle versioneer
from ._version import get_versions
__author__ = """Paul Saxe"""
__email__ = 'psaxe@molssi.org'
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

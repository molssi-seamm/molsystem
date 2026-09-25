***************
Getting Started
***************

Installation
============
MolSystem is probably already installed in your SEAMM environment, but
if not or if you wish to check, follow the directions for the `SEAMM Installer`_. The
graphical installer is the easiest to use. In the SEAMM conda environment, simply type::

  seamm-installer

or use the shortcut if you installed one. Switch to the second tab, `Components`, and
check for `quickmin-step`. If it is not installed, or can be updated, check the box
next to it and click `Install selected` or `Update selected` as appropriate.

The non-graphical installer is also straightforward::

  seamm-installer install --update molsystem

will ensure both that it is installed and up-to-date.

MolSystem is not directly tied to SEAMM, so can be used outside SEAMM. Install it
from PyPI with pip (or uv) into any Python 3.10 or later environment:

  pip install molsystem

Its compiled dependencies, RDKit and Open Babel, are installed automatically from
their PyPI wheels on macOS, Linux and Windows. A conda-forge package also exists:

  conda install -c conda-forge molsystem

but it lags the PyPI release, and you should not mix the two in one environment.

.. _SEAMM Installer: https://molssi-seamm.github.io/installation/index.html

Replace this!
=============
Put an example or two here....

That should be enough to get started. For more detail about the functionality in this plug-in, see the :ref:`User Guide <user-guide>`.

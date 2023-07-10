=============
MolSystem
=============

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/dftbplus_step
   :target: https://github.com/molssi-seamm/molsystem/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/molsystem/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/molsystem/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/molsystem/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/molsystem
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/molsystem/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/molsystem/security/code-scanning
   :alt: Code Quality


.. image:: https://github.com/molssi-seamm/molsystem/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/molsystem/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/molsystem.svg
   :target: https://pypi.python.org/pypi/molsystem
   :alt: PyPi VERSION

Description
-----------
Molsystem provides a general class for handling molecular and periodic systems. This is
the heart of SEAMM.

* Free software: GNU Lesser General Public License v3+
* Documentation: https://molsystem.readthedocs.io.


Features
--------

* Supports molecular and periodic (crystalline) systems.
* Support for pointgroup and spacegroup symmetry.
* Implemented as a SQL database, currently using SQLite, which provides permanence and a
  disk file.
* Handles multiple systems, and for each system any number of
  configurations. Configurations are different structures, generally having different
  coordinates. Examples are conformers or frames in an MD trajectory.
* Configurations can have differing bonds, supporting reactive forcefields such as
  ReaxFF.
* Configurations can also have differing numbers of atoms, supporting e.g. grand
  canonical Monte Carlo.
* Support for templates and subsets. Templates can be used for creating or finding
  structures in systems and subsets hold substes of atoms. An atom can be in multiple
  subsets. For proteins, subsets can be used for residues and chains; for polymers,
  monomers; and for fluids the subsets can track individual molecules.
* Subsets and templates handle the correspondance between atoms so if, for instance, the
  ordering of atoms in residues are different, that can be handled using the
  correspondances.
* Direct connection to and support of RDKit and OpenBabel molecules, allowing direct
  transfer back and forth between configurations and those libraries.
* Properties of systems and configurations are stored in fact tables of integers,
  floats, strings, and JSON, using a star schema. One of the dimensions is metadata
  describing the property. This data warehousing approach is almost identical to that
  used in CIF files, with the metadata being the CIF dictionaries.
* Good support of CIF files, with the ability to store data other than the actual atoms
  using the properties warehouse.
* Direct support for SMILES, SMARTS and substructure searching.
* Scales to millions of systems and configurations, and databases in excess of 100 GB.

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award OAC-1547580 and CHE-2136142.

.. _MolSSI: https://www.molssi.org
.. _`National Science Foundation`: https://www.nsf.gov

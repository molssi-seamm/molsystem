=======
History
=======
2023.10.30 -- Support for InChI improved, RMSD and PubChem added...
    * Adds support for aligning structures and calculating RMSD
    * Adds support for working directly with PubChem to get structures, IUPAC names,
      etc.
    * Improves support for InChI, working around issues in both OpenBabel and RDKit.
    * Added substantial new functionality for spacegroups and primitive cell handling,
      but still not complete.

2023.9.20 -- Better support for primitive cells and spacegroups
    * Added getting the spacegroup from the symmetry operators
    * Fixed updating the coordinates from the primitive cell

2023.9.5 -- Support for velocities of atoms.

2023.8.30 -- Support for spacegroup symmetry.

2023.8.27 -- Bugfix: writing SDF did not handle charge and multiplicity.

2023.7.30 -- Improved handling of properties
    * Added ability to get lists of systems or configurations filtered by name
    * Improved handling of properties on just a system, not configuration
    * Added ability to filter properties retrieved
    * Improved handling of properties when creating OpenBabel OB_MOL object
      
2023.7.26 -- Bugfix: error in QCSchema bonds; enhancement: RDKit
    * Fixed bug in the bond indices in QCSchema
    * Added ability to use RDKit for SMILES and InChI

2023.7.18.1 -- Added support for creating structures from InChIKeys
    * Uses PubChem to translate the InChiKey to InChI.
       
2023.7.18 -- Added support for InChI and InChIKeys

2023.7.9 -- Added JSON properties
    * Added properties stored as JSON, which allows, vectors, tensors, etc.
      
2023.4.6 -- Enhancements for CIF files
    * Handle uncertainties in CIF files expressed as '(x)' at end of value.

2023.3.30 -- Enhancements to QCSchema support
    * Improved naming of molecule in QCSchema
    * Added creation of configurations from QCSchema objects.

2023.2.13 -- Fixed issue with valence in RDkit for cations like NH4+

2022.11.20 -- Added a method to copy a configuration.
  Added a new method to the `system` class, `copy_configuration`, that creates a copy of
  the configuration using the same atomset and bonset, but new coordinates and cell so
  that any changes to coordinates and cell are not shared between the configurations. By
  default it copies the current configuration.

2022.11.18 -- Fixed bug with handling for Open Babel
  The total charge and multiplicity were not correctly set when creating an Open Babel
  molecule.

2022.11.3 -- Add handling of strain and improved handling of properties
  Added methods for straining the unit cell, and also straining a configuration,
  correctly handling the coordinates for an affine transformation. In the future will
  add e.g. affine transformation of the centers of molecules, which is useful for
  molecular fluids.

  Added the system for properties, in addition to the configuration. This allows system
  properties that are not associated with a particular configuration, which is often
  appropriate for experimental results. It also makes it much easier to search for
  systems where any configuration has a particular property.

2022.10.26 -- Improved database write performance.
  Switched to write-ahead mode and tweaked memory settings. This gives a large
  performance improvement (10x or more) for large database (~1 GB).

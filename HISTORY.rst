=======
History
=======

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

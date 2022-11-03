=======
History
=======

2022.11.3 -- Add handling of strain.
  Added methods for straining the unit cell, and allso straining a configuration,
  correctly handling the coordinates for an affine transformation. In the future will
  add e.g. affine transformation of the centers of molecules, which is useful for
  molecular fluids.

2022.10.26 -- Improved database write performance.
  Switched to write-ahead mode and tweaked memory settings. This gives a large
  performance improvement (10x or more) for large database (~1 GB).

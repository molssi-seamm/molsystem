.. _user-guide:

**********
User Guide
**********
MolSystem is a Python package for storing molecular and periodic (crystalline)
systems. It provides and object oriented interface over an underly SQL database that
handles systems with multiple configurations, plus templates and subsets. You can
loosely think of a system as being an abstract molecular or periodic system and the
configurations as the physical instances with the coordinates, bonds, etc. that describe
different structures, conformers, or frames in a trajectory.

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Perceiving bonds
================
Structures read from formats without connectivity -- an extended XYZ trajectory frame,
for example -- have atoms but no bonds, so anything that works with molecules (finding
molecules, keeping them whole across a periodic boundary, extracting clusters) has
nothing to go on. ``perceive_bonds()`` finds the bonds from the geometry and adds them to
the configuration::

    n_bonds = configuration.perceive_bonds()

Two atoms are bonded when their distance is less than ``tolerance`` (default 1.2) times
the sum of their covalent radii (Pyykkö radii from ``mendeleev``). Hydrogen is limited
to one bond, the shortest, so a hydrogen bond is never mistaken for a covalent one, and
the alkali and alkaline-earth metals are treated as ions and left unbonded. Periodic
configurations use the minimum image, so molecules straddling the cell boundary are
bonded correctly. All bonds are single bonds; bond orders are not assigned.

The defaults can be adjusted per call::

    configuration.perceive_bonds(
        replace=True,            # discard existing bonds first
        tolerance=1.25,          # looser criterion, e.g. hot MD frames
        radii={"H": 0.35},       # override a covalent radius (Å)
        max_bonds={"H": 1, "O": 2},  # per-element valence limits
        exclude=(),              # bond the metals too, e.g. an ionic crystal
    )

Only the ``covalent radii`` method is available at present; the ``method`` argument is
there so that others, such as a Voronoi tessellation, can be added later.

Index
=====

* :ref:`genindex`

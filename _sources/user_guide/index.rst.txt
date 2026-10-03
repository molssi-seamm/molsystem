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
the alkali and alkaline-earth metals are treated as ions and left unbonded. In a
periodic configuration every periodic image within reach is considered and each bond
records the cell offset of its partner, so molecules straddling the cell boundary are
bonded correctly and covalent crystals get all their bonds -- primitive diamond, for
example, has four bonds between its two atoms through different images. All bonds are
single bonds; bond orders are not assigned.

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

Charges and spin in structure files
===================================
A configuration carries its charge and spin multiplicity as properties of the structure
as a whole::

    configuration.charge = -1
    configuration.spin_multiplicity = 1

The atoms may *also* carry formal charges of their own, in the optional
``formal_charge`` attribute -- the charge on the carboxylate oxygen of an acetate ion,
for instance, rather than on the ion as a whole::

    configuration.atoms.get_column_data("formal_charge")   # e.g. [0, 0, 0, -1]

The two are distinct and are written to a structure file as such: the molecular charge
and multiplicity as properties of the structure, and the formal charges on the
individual atoms that carry them -- the ``M  CHG`` line of an SDF file, say. Both come
back unchanged on reading the file, so the charges stay where they belong and a neutral,
closed-shell structure is not mistaken for a radical.

Properties and their units
==========================
Every property has units, held as a string that Pint can interpret, and a property that
is genuinely dimensionless -- a statistical inefficiency, a count, a ratio -- has an
*empty* unit string rather than no units at all::

    configuration.properties.units("temperature, inefficiency#LAMMPS#oplsaa+")   # ""

The distinction matters because the units are what a value is converted to when a
property is stored or read, so a property whose units are undefined cannot be converted.
Anything that creates properties -- including reading them back from a structure file --
records dimensionless units as ``""``.

User tables
===========
Besides the structures, the database holds the tables a SEAMM flowchart builds -- with
the Table step, or by storing results into a table -- so that a job's tables are in its
``seamm.db`` with everything else. They are reached by name through
``SystemDB.user_tables``::

    tables = db.user_tables
    table = tables.create(
        "energies",
        columns=[("SMILES", "string", None), ("E (kJ/mol)", "float", None)],
        index_column="SMILES",
    )
    row = table.append_row(SMILES="CCO")
    table.set_cell(row, "E (kJ/mol)", -235.1)
    df = table.to_dataframe()          # a pandas copy, typed by the declared types

Each column has a declared type -- ``boolean``, ``integer``, ``float``, ``string`` or
``json`` -- and a default, which fills the column in rows that do not give a value. The
values are stored exactly as written (the text ``"1.0960"`` stays text) and read back
according to the declared type; values written to text columns are stored as text. A
table may have an *index column* whose values identify the rows; rows also have an
internal id, which is not meant for users, and a *current row* that SEAMM's steps write
to (``None`` means the next write appends a row).

The tables are stored under internal names, ``table_1``, ``table_2``, ..., with
internal column names, so any display name is allowed, including names that differ only
in case. A registry table, ``_tables``, records each table's name, columns, index column,
current row and metadata, and a journal, ``_table_changes``, records every change
(creating and dropping tables, adding columns, appending rows and setting values). Row
ids are never reused. Nothing is committed by these calls; the caller commits.

A database opened read-only (``file:seamm.db?mode=ro``) can be read but not written:
creating or changing a table raises ``PermissionError``.


Index
=====

* :ref:`genindex`

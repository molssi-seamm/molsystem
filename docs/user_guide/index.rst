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
single bonds: anything that derives a SMILES or Lewis structure from them would read
ethylene as ``[CH2+][CH2-]``. For a molecule, ``assign_bond_orders()`` then gives the
bonds their orders and the atoms their formal charges::

    n_assigned, n_failed = configuration.assign_bond_orders(charge=0)

Each molecule (connected fragment) without metals gets the bond orders of RDKit's
``DetermineBondOrders``. The first closed-shell structure is taken, trying the molecule's
own charge when it is the only molecule, and otherwise the charges 0, -1, +1, -2 and +2
in turn; a structure without formal charges is tried before one with them (N#As, not
[N+]#[As-]). Failing that, a structure with as many radical electrons as the spin
multiplicity says is taken (the ethyl radical; the vinyl radical, from the closed-shell
anion). A molecule with neither keeps its single bonds, with radical electrons where its
atoms lack valence: Ge2 is [Ge][Ge] and GeC [C][Ge], not charge-separated. A molecule
RDKit cannot handle at all, such as hypervalent PF6-, keeps its single bonds, which are
already its Lewis structure. Single alkali, alkaline-earth and halogen atoms get their
ionic charges. Given the net charge, the molecules' charges are checked against it and
a mismatch is logged. The Lewis structure does not depend on the spin multiplicity:
triplet AsN is N#[As]. Periodic configurations are left alone. Read Structure does this
for extended XYZ and MOPAC files when it perceives the bonds.

The radical electrons are stored with the atoms, in the optional ``radical_electrons``
attribute (as the formal charges are in ``formal_charge``), and SMILES and RDKit
molecules use them as they are rather than fitting the charges and radicals to the spin
multiplicity. SDF files hold them as ``M  RAD`` for one or two radical electrons on an
atom and, since MDL files cannot hold more (three on each Ge of [Ge][Ge]), also as the
data item ``SEAMM|radical electrons|json|``.

Metal complexes
---------------
Bonds to metals (``molsystem.topology.METALS``: the alkali, alkaline-earth and
transition metals, lanthanides, actinides and the heavier post-transition metals) follow
one convention, ionic + dative:

* every metal-ligand bond is a dative bond from the ligand's donor atom to the metal,
  bond order 9 (MDL's bond type 9, RDKit's ``DATIVE``);
* each ligand -- a fragment of the complex without its metals -- has its closed-shell
  Lewis structure with the least charge separation: CO neutral (``[C-]#[O+]``), Cp-,
  Cl-, acetate, NMe2-; a lone atom is the anion of its default valence: H-, O2-, N3-,
  C4-, Si4-; a pi ligand (Cp, arene, diene) has a dative bond from each of its bonded
  atoms;
* metal-metal bonds are single bonds;
* the metals carry the rest of the net charge: the alkali metals +1 and alkaline-earth
  metals +2, the others sharing what is left as evenly as integers allow. Fe(CO)5 is
  Fe(0), ferrocene Fe2+, Cp2TiCl2 Ti4+, ZrC Zr4+ and C4-, NaZrF5 Na+ and Zr4+. When
  only alkali and alkaline-earth metals are present, the ligands take the charges that
  balance them if they can (Cs2O2: the peroxide O2 2-, not O=O); otherwise these metals
  give way (SrOH: Sr+). The same holds for their ions and the molecules beside them,
  since perception leaves these metals unbonded.

The metals' formal charges are bookkeeping, not always oxidation states. The unpaired
electrons of an open-shell complex are the metal's, so complexes are assigned whatever
their spin multiplicity, and their ligands carry no radicals.

``perceive_bonds()`` uses a wider tolerance for pairs with a metal (``metal_tolerance``,
default 1.3), since bonds to pi ligands are longer than the covalent radii suggest. It
drops a metal's bond to the central atom of a four-membered chelate ring -- the C of a
kappa-2 carboxylate, the N of a nitrate -- which that tolerance reaches, when it is at
least 0.3 Å longer than the bonds to both chelating atoms. Distance-based perception
can still be wrong for long metal-ligand contacts, agostic hydrogens, bridging hydrides,
ion pairs and clusters: check such structures.

Dative bonds keep their direction. Bonds are stored with the first atom's index the
lower, so a dative bond whose donor is the second atom has bond order 10
(``molsystem.DATIVE_REVERSED``); ``molsystem.dative_donor_acceptor(i, j, bondorder)``
gives (donor, acceptor). SDF and MDL files, Open Babel and RDKit all get the donor
first.

The defaults can be adjusted per call::

    configuration.perceive_bonds(
        replace=True,            # discard existing bonds first
        tolerance=1.25,          # looser criterion, e.g. hot MD frames
        metal_tolerance=1.3,     # for pairs with a metal
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

Radicals and ions in SMILES and RDKit
-------------------------------------
``to_smiles()`` and ``to_RDKMol()`` carry the charges and unpaired electrons into RDKit,
so that radicals and ions give the right SMILES::

    configuration.from_smiles("[CH2]C")        # charge 0, multiplicity 2
    configuration.to_smiles(canonical=True)    # '[CH2]C', the ethyl radical

``from_smiles()`` sets the charge and spin multiplicity from the SMILES, and the formal
charges of the atoms. Since the hydrogens are all explicit, an atom short of its valence
is a radical, an anion or a cation: the formal charges are used when they add up to the
configuration's charge, and otherwise, as for a structure read from a file without them,
charges and radical electrons are chosen that match the charge and multiplicity, with
the fewest charges, negative charges on electronegative atoms, positive ones on metals
and electropositive atoms, and charge-separated groups such as the nitro group
C[N+](=O)[O-] kept neutral. Hypervalent atoms need no charge of their own: sulfuric acid
is O=S(=O)(O)O. The same atoms and bonds are the radical [CH2]C, the cation [CH2+]C or the
anion [CH2-]C, depending on the configuration's charge and multiplicity.

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

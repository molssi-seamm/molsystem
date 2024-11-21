#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fixtures for testing the 'molsystem' package."""
import math
from pathlib import Path

import pytest

from molsystem.system_db import SystemDB

path = Path(__file__).resolve().parent
data_path = path / "data"


def pytest_addoption(parser):
    parser.addoption(
        "--run-timing", action="store_true", default=False, help="run timing tests"
    )
    parser.addoption(
        "--openeye", action="store_true", default=False, help="run openeye tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "timing: mark test as timing to run")
    config.addinivalue_line(
        "markers", "openeye: mark test as using OpenEye functionality"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-timing"):
        # --run-timing given in cli: do not skip timing tests
        return
    skip_timing = pytest.mark.skip(reason="need --run-timing option to run")
    for item in items:
        if "timing" in item.keywords:
            item.add_marker(skip_timing)

    if config.getoption("--openeye"):
        # --openeye given in cli: include the openeye tests
        return
    skip_openeye = pytest.mark.skip(reason="need --openeye option to run")
    for item in items:
        if "openeye" in item.keywords:
            item.add_marker(skip_openeye)


def mk_table(db, name="table1"):
    """Create a table with some attributes."""
    table = db.create_table(name)
    table.add_attribute("atno", coltype="int", default=-1)
    for column in ("x", "y", "z"):
        table.add_attribute(column, coltype="float")
    return table


@pytest.fixture()
def empty_db():
    """Create a system db with no systems."""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")

    yield db

    db.close()
    try:
        del db
    except:  # noqa: E722
        print("Caught error deleting the database")


@pytest.fixture()
def db(empty_db):
    """Create an empty system db."""
    system = empty_db.create_system(name="default")
    system.create_configuration(name="default")

    return empty_db


@pytest.fixture()
def system(db):
    """An empty system."""
    return db.system


@pytest.fixture()
def configuration(system):
    """An empty system."""
    return system.configuration


@pytest.fixture()
def symmetry(configuration):
    """An _Symmetry object."""
    return configuration.symmetry


@pytest.fixture()
def atoms(configuration):
    """An empty atoms table."""
    return configuration.atoms


@pytest.fixture()
def bonds(configuration):
    """An empty bonds table."""
    return configuration.bonds


@pytest.fixture()
def two_dbs():
    """Create two different dbs."""
    db1 = SystemDB(filename="file:seamm_db1?mode=memory&cache=shared")
    system1 = db1.create_system(name="default")
    system1.create_configuration(name="default")

    db2 = SystemDB(filename="file:seamm_db2?mode=memory&cache=shared")
    system2 = db2.create_system(name="default")
    system2.create_configuration(name="default")

    yield db1, db2

    db1.close()
    db2.close()


@pytest.fixture()
def two_systems(two_dbs):
    """Two empty systems."""
    db1, db2 = two_dbs
    return db1.system, db2.system


@pytest.fixture()
def two_configurations(two_systems):
    """Two empty configuration."""
    system1, system2 = two_systems
    return system1.configuration, system2.configuration


@pytest.fixture()
def simple_table(db):
    return mk_table(db)


@pytest.fixture()
def two_tables(db):
    return mk_table(db, "table1"), mk_table(db, "table2")


@pytest.fixture()
def system_with_two_tables(db):
    mk_table(db, "table1")
    mk_table(db, "table2")
    return db.system


@pytest.fixture()
def db_with_two_tables(db):
    mk_table(db, "table1")
    mk_table(db, "table2")
    return db


@pytest.fixture()
def two_tables_in_two_dbs(two_dbs):
    db1, db2 = two_dbs
    return mk_table(db1, "table1"), mk_table(db2, "table2")


@pytest.fixture()
def old_templates(system):
    """A system with a template for water."""
    templates = system["template"]
    atoms = system["templateatom"]
    atoms.add_attribute("name", "str", default="")
    bonds = system["templatebond"]

    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    name = ["O", "H1", "H2"]
    i_atom = [0, 0]
    j_atom = [1, 2]

    tid = templates.append(name="H2O", type="molecule")[0]
    templates.current_template = tid

    ids = atoms.append(atno=atno, x=X, y=Y, z=Z, name=name)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    bonds.append(i=j, j=i)  # flipped on purpose so code orders.

    # Acetic acid

    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    X = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    Y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    Z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    symbol = ['C', 'H', 'H', 'H', 'C', 'O', 'O', 'H']
    name = ['C1', 'H1', 'H2', 'H3', 'C', 'O', 'OH', 'H']

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    tid = templates.append(name="acetic acid", type="molecule")[0]
    templates.current_template = tid

    ids = atoms.append(symbol=symbol, x=X, y=Y, z=Z, name=name)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    bonds.append(i=i, j=j, bondorder=order)

    return system


@pytest.fixture
def H2O(configuration):
    """A configuration of  water."""
    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    name = ["O", "H1", "H2"]
    i_atom = [0, 0]
    j_atom = [1, 2]

    configuration.system.name = "water"
    configuration.name = "TIP3P"
    configuration.atoms.add_attribute("name", "str", default="")

    ids = configuration.atoms.append(x=X, y=Y, z=Z, atno=atno, name=name)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    configuration.bonds.append(i=i, j=j, bondorder=1)

    return configuration


@pytest.fixture()
def AceticAcid(configuration):
    """An configuration object for an acetic acid molecule"""
    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    x = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    atno = [6, 1, 1, 1, 6, 8, 8, 1]  # noqa: E221

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    configuration.system.name = "acetic acid"
    configuration.name = "acetic acid"
    ids = configuration.atoms.append(x=x, y=y, z=z, atno=atno)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    configuration.bonds.append(i=i, j=j, bondorder=order)

    return configuration


@pytest.fixture()
def Acetate(configuration):
    """An configuration object for an acetate acid molecule"""
    # yapf: disable
    #       C       H        H        H        C        =O      O
    x = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757]  # noqa: E221, E501, E201
    y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970]  # noqa: E221, E501, E201
    z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919]  # noqa: E221, E501
    atno = [6, 1, 1, 1, 6, 8, 8]  # noqa: E221

    #       C  H  H  H  C =O  O
    i_atom = [0, 0, 0, 0, 4, 4]
    j_atom = [1, 2, 3, 4, 5, 6]
    order =  [1, 1, 1, 1, 2, 1]  # noqa: E222
    # yapf: enable

    configuration.system.name = "acetate"
    configuration.name = "acetate"
    ids = configuration.atoms.append(x=x, y=y, z=z, atno=atno)

    i = [ids[x] for x in i_atom]
    j = [ids[x] for x in j_atom]

    configuration.bonds.append(i=i, j=j, bondorder=order)
    configuration.charge = -1

    properties = configuration.properties
    properties.add("float property", "float", units="kcal/mol")
    properties.add("int property", "int")
    properties.add("str property", "str")

    properties.put("float property", 3.14)
    properties.put("int property", 2)
    properties.put("str property", "Hi!")

    return configuration


@pytest.fixture
def disordered(AceticAcid):
    """Two acetic acid molecules with different atom order"""
    # yapf: disable
    #       C       H        H        H        C        =O      O        H
    x = [ 1.0797, 0.5782,  0.7209,  0.7052,  0.5713, -0.1323, 0.9757,  2.1724]  # noqa: E221, E501, E201
    y = [ 0.0181, 3.1376, -0.6736, -0.3143,  1.3899,  1.7142, 2.2970,  0.0161]  # noqa: E221, E501, E201
    z = [-0.0184, 0.2813, -0.7859,  0.9529, -0.3161, -1.2568, 0.5919, -0.0306]  # noqa: E221, E501
    atno = [6, 1, 1, 1, 6, 8, 8, 1]  # noqa: E221

    #       C  H  H  H  C =O  O  H
    i_atom = [0, 0, 0, 0, 4, 4, 6]
    j_atom = [1, 2, 3, 4, 5, 6, 7]
    order =  [1, 1, 1, 1, 2, 1, 1]  # noqa: E222
    # yapf: enable

    configuration = AceticAcid
    ids = configuration.atoms.append(
        x=[*reversed(x)], y=[*reversed(y)], z=[*reversed(z)], atno=[*reversed(atno)]
    )

    i = [ids[7 - x] for x in i_atom]
    j = [ids[7 - x] for x in j_atom]

    configuration.bonds.append(i=i, j=j, bondorder=order)

    return configuration


@pytest.fixture()
def vanadium(configuration):
    """BCC vanadium crystal, without symmetry."""
    configuration.name = "BCC Vanadium"
    configuration.periodicity = 3
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = [3.03, 3.03, 3.03, 90, 90, 90]
    configuration.atoms.append(
        x=[0.0, 0.5],
        y=[0.0, 0.5],
        z=[0.0, 0.5],
        gx=[0.0, 0.02],
        gy=[0.01, -0.02],
        gz=[-0.01, 0.01],
        vx=[0.0, 0.02],
        vy=[0.01, -0.02],
        vz=[-0.01, 0.01],
        symbol="V",
    )
    return configuration


@pytest.fixture()
def copper(configuration):
    """FCC copper crystal, without symmetry."""
    x = [0.0, 0.5, 0.5, 0.0]
    y = [0.0, 0.5, 0.0, 0.5]
    z = [0.0, 0.0, 0.5, 0.5]
    configuration.name = "FCC Copper"
    configuration.periodicity = 3
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (3.61491, 3.61491, 3.61491, 90, 90, 90)
    configuration.atoms.append(x=x, y=y, z=z, symbol=["Cu"])
    return configuration


@pytest.fixture()
def diamond(configuration):
    """Diamond crystal."""
    # Atoms
    x = [0.0]
    y = [0.0]
    z = [0.0]

    # Bonds
    # Generators = [ 0,  2,  1,  5, 10, 18,  3,  7]
    II = [0, 0, 0, 0]
    JJ = [0, 0, 0, 0]
    symop2 = ["2", "6_544", "4_454", "8_445"]

    configuration.name = "Diamond"
    configuration.periodicity = 3
    configuration.group = "F d -3 m :1"
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (3.57, 3.57, 3.57, 90, 90, 90)
    ids = configuration.atoms.append(x=x, y=y, z=z, symbol=["C"])

    Is = [ids[i] for i in II]
    Js = [ids[j] for j in JJ]
    configuration.bonds.append(i=Is, j=Js, symop2=symop2)

    return configuration


@pytest.fixture()
def h_chain(configuration):
    """Hydrogen chain crystal."""
    # Atoms
    x = [0.0]
    y = [0.0]
    z = [0.5]

    # Bonds
    II = [0]
    JJ = [0]
    symop2 = ["1_556"]

    configuration.name = "H-chain"
    configuration.periodicity = 3
    configuration.group = "P 1"
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (10.0, 10.0, 1.0, 90, 90, 90)
    ids = configuration.atoms.append(x=x, y=y, z=z, symbol=["H"])

    Is = [ids[i] for i in II]
    Js = [ids[j] for j in JJ]
    configuration.bonds.append(i=Is, j=Js, symop2=symop2)

    return configuration


@pytest.fixture()
def h_chain2(configuration):
    """Hydrogen chain crystal."""
    # Atoms
    x = [0.0, 0.0]
    y = [0.0, 0.0]
    z = [0.26, 0.74]

    # Bonds
    II = [0, 0]
    JJ = [1, 1]
    symop2 = [".", "1_554"]

    configuration.name = "H-chain2"
    configuration.periodicity = 3
    configuration.group = "P 1"
    configuration.coordinate_system = "fractional"
    configuration.cell.parameters = (10.0, 10.0, 2.0, 90, 90, 90)
    ids = configuration.atoms.append(x=x, y=y, z=z, symbol=["H", "H"])

    Is = [ids[i] for i in II]
    Js = [ids[j] for j in JJ]
    configuration.bonds.append(i=Is, j=Js, symop2=symop2)

    return configuration


@pytest.fixture()
def CH3COOH_3H2O(AceticAcid):
    """Configuration with acetic acid and 3 water molecules"""
    configuration = AceticAcid

    # TIP3P
    r0 = 0.9572
    theta0 = 104.52

    # H locations are ±x, 0, z
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))

    X = [0.0, x, -x]
    Z = [0.0, z, z]

    atno = [8, 1, 1]
    i_atom = [0, 0]
    j_atom = [1, 2]

    configuration.name = "acetic acid with 3 waters"

    for no in range(1, 4):
        ids = configuration.atoms.append(x=X, y=no * 5.0, z=Z, atno=atno)

        i = [ids[x] for x in i_atom]
        j = [ids[x] for x in j_atom]

        configuration.bonds.append(i=i, j=j)

    return configuration


# Fixtures for working with the files in data/


@pytest.fixture()
def full_db():
    """Create a system database with several systems."""
    db = SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    for filename in ["acy.mmcif"]:
        with open(data_path / filename, "r") as fd:
            text = fd.read()
        system = db.add_system(name=Path(filename).stem)
        configuration = system.create_configuration(name="default")
        configuration.from_mmcif_text(text)

    yield db

    db.close()


@pytest.fixture(scope="session")
def amino_acids():
    """Create a system database with 20 amino acids, each as a system."""
    db = SystemDB(filename="file:amino_acids_db?mode=memory&cache=shared")
    db.read_cif_file(data_path / "aminoacids.mmcif")

    yield db

    db.close()


@pytest.fixture(scope="session")
def aa_templates():
    """Create a system database with 20 amino acids, each as a template."""
    db = SystemDB(filename="file:aa_templates_db?mode=memory&cache=shared")
    db.read_cif_file(data_path / "aminoacids.mmcif")
    templates = db.templates
    for system in db.systems:
        templates.create(
            name=system.name,
            category="amino acid",
            configuration=system.configuration.id,
        )

    yield templates

    db.close()


@pytest.fixture()
def gly(aa_templates):
    """The template for glycine."""
    return aa_templates.get(8)


@pytest.fixture()
def simple_templates(CH3COOH_3H2O):
    """Simple templates for water and acetic acid, plus solvated acetic acid."""
    db = CH3COOH_3H2O.system_db
    templates = db.templates

    db.read_cif_file(data_path / "acy.mmcif")[0]
    db.read_cif_file(data_path / "hoh.mmcif")[0]

    templates.create(name="acy", category="molecule")
    templates.create(name="hoh", category="molecule")

    return db


@pytest.fixture()
def full_templates(CH3COOH_3H2O):
    """Full templates for water and acetic acid, plus solvated acetic acid."""
    db = CH3COOH_3H2O.system_db
    templates = db.templates
    system = CH3COOH_3H2O.system
    system.name = "acetic acid"

    acy_sys = db.read_cif_file(data_path / "acy.mmcif")[0]
    acy_conf = db.get_system(acy_sys).configuration
    acy_cid = acy_conf.id
    hoh_sys = db.read_cif_file(data_path / "hoh.mmcif")[0]
    hoh_conf = db.get_system(hoh_sys).configuration
    hoh_cid = hoh_conf.id

    templates.create(name="acy", category="molecule", configuration=acy_cid)
    templates.create(name="hoh", category="molecule", configuration=hoh_cid)

    return db


@pytest.fixture()
def properties(system):
    """Lots of configurations with many properties."""
    properties = system.system_db.properties

    # Create systems
    cids = []
    for i in range(1000):
        configuration = system.create_configuration(name=str(i))
        cids.append(configuration.id)

    # And properties
    for i in range(10):
        fid = properties.add(f"float_{i}", "float")
        iid = properties.add(f"int_{i}", "int")
        sid = properties.add(f"str_{i}", "str")

        fval = -10.0 + i
        ival = -500 + 100 * i
        for cid in cids:
            properties.put(cid, fid, fval)
            properties.put(cid, iid, ival)
            properties.put(cid, sid, f"string {ival}")
            fval += 0.01
            ival += 1

    return properties


@pytest.fixture()
def polyethylene(configuration):
    """A polyethylene crystal with bonds."""
    path = data_path / "polyethylene.cif"
    cif_text = path.read_text()
    configuration.from_cif_text(cif_text)

    return configuration


@pytest.fixture()
def benzene(configuration):
    """A benzene crystal with bonds."""
    path = data_path / "benzene.cif"
    cif_text = path.read_text()
    configuration.from_cif_text(cif_text)

    return configuration

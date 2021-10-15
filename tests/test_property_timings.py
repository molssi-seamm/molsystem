#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for properties."""

from pathlib import Path  # noqa: F401
import pprint  # noqa: F401
import pytest  # noqa: F401
import time

from molsystem import SystemDB  # noqa: F401

nsystems = 1000000
nproperties = 10


@pytest.fixture(scope="module")
def db_file(tmp_path_factory):
    """Create and temporary directory, and clean it up at the end."""
    dirpath = tmp_path_factory.mktemp("data")
    path = dirpath / "seamm.db"
    db = SystemDB(filename=str(path))
    db.close()

    print(f"temp file: {path}")

    yield str(path)

    path.unlink()
    dirpath.rmdir()


@pytest.fixture()
def disk_db(db_file):
    """Create and return a SystemDB object with the database file."""
    db = SystemDB(filename=db_file)

    yield db

    db.close()


@pytest.fixture()
def saved_db():
    """Create and return a SystemDB object with the database file."""
    db = SystemDB(filename="seamm.db")

    yield db

    db.close()


@pytest.fixture()
def properties(disk_db):
    """Lots of configurations with many properties."""
    # path = Path("seamm.db")
    # exists = path.exists()

    # db = SystemDB(filename="seamm.db")
    exists = False
    db = disk_db

    float_data = db["float_data"]
    int_data = db["int_data"]
    str_data = db["str_data"]

    properties = db.properties

    if not exists:
        t0 = time.perf_counter()
        # Create systems
        cids = []
        names = [str(i) for i in range(nsystems)]
        systems = db["system"].append(name=names)
        cids = db["configuration"].append(system=systems, name="1")

        t1 = time.perf_counter()
        print(f"  {nsystems} systems took {t1-t0:.3f} s")

        # And properties
        fstart = -nproperties
        delta = 10.0 / nsystems
        istart = -nsystems
        for i in range(nproperties):
            fid = properties.add(f"float_{i}", "float")
            iid = properties.add(f"int_{i}", "int")
            sid = properties.add(f"str_{i}", "str")

            fval = fstart + i
            ival = istart + nsystems * i
            fvals = []
            ivals = []
            svals = []
            for cid in cids:
                fvals.append(fval)
                ivals.append(ival)
                svals.append(f"string {ival}")
                fval += delta
                ival += 1
            float_data.append(configuration=cids, property=fid, value=fvals)
            int_data.append(configuration=cids, property=iid, value=ivals)
            str_data.append(configuration=cids, property=sid, value=svals)

        t1 = time.perf_counter()
        print(f"  db creation took {t1-t0:.3f} s")

    yield properties

    db.close()


@pytest.mark.timing
def test_query2(properties):
    """Test a simple query"""

    t0 = time.perf_counter()
    result = properties.query(
        "float_0",
        "between",
        -10.0,
        -10.0 + 1000 / nsystems,
        what=["float_8", "str_9"],
    )
    t1 = time.perf_counter()
    print(f"  test took {t1-t0:.6f} s")
    print(f"        retrieved {len(result['str_9'])} strings")
    print(f"              and {len(result['float_8'])} floats")

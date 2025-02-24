#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `molsystem` package.

What is the version of SQLITE?
"""
import sqlite3


def test_sqlite_version():
    """Print the version of SQLite"""
    print(f"SQLite version {sqlite3.sqlite_version}")
    assert True

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

import molsystem
"""Tests for the system class."""


def test_construction():
    """Simplest test that we can make a System object"""
    system = molsystem.System()
    assert str(type(system)) == "<class 'molsystem.system.System'>"


def test_version_empty():
    """Simplest test that we can make a System object"""
    system = molsystem.System()
    assert system.version == 0


def test_version_changed():
    """Simplest test that we can make a System object"""
    system = molsystem.System()
    with system as sys:
        sys.periodicity = 3
    assert system.version == 1


def test_version_unchanged():
    """Simplest test that we can make a System object"""
    system = molsystem.System()
    with system as sys:  # noqa: F841
        pass
    assert system.version == 0

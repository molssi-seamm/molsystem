#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `bonds` in the `molsystem` package."""

import pytest  # noqa: F401


def test_simple(AceticAcid):
    """Simple test for acetic acid"""
    assert AceticAcid.formula == ("C2 H4 O2", "C H2 O", 2)


def test_periodic(copper):
    """Test for copper crystal"""
    assert copper.formula == ("Cu4", "Cu", 4)


def test_charged(configuration):
    """Test for charged configuration"""
    configuration.from_smiles("[NH4+]")
    assert configuration.formula == ("[H4 N]+", "[H4 N]+", 1)


def test_doubly_charged(configuration):
    """Test for charged configuration"""
    configuration.from_smiles("[NH4+].[NH4+]")
    assert configuration.formula == ("[H8 N2]+2", "[H4 N]+", 2)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

"""Tests for handling PubChem."""


def test_to_smiles(AceticAcid):
    """Get the name from a system"""
    correct = "acetic acid"
    name = AceticAcid.PC_iupac_name()

    if name != correct:
        print(name)
    assert name == correct


def test_from_cid(configuration):
    """Get the system from a PubChem cid"""
    correct = "(Z)-but-2-ene"
    cid = 5287573
    configuration.PC_from_cid(cid)

    new_cid = configuration.PC_cid
    name = configuration.PC_iupac_name()

    if new_cid != cid:
        print(new_cid)
    assert new_cid == cid

    if name != correct:
        print(name)
    assert name == correct


def test_from_name(configuration):
    """Get the system from a name"""
    correct = "(Z)-but-2-ene"
    cid = 5287573
    configuration.PC_from_identifier("cis-2-butene")

    new_cid = configuration.PC_cid
    name = configuration.PC_iupac_name()

    if new_cid != cid:
        print(new_cid)
    assert new_cid == cid

    if name != correct:
        print(name)
    assert name == correct

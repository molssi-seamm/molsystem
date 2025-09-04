#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the QCSchema mixin of the class."""

import pprint  # noqa: F401
import pytest  # noqa: F401


def test_water(H2O):
    """Test creating QCSchema for water."""
    correct = '{"schema_name": "qcschema_molecule", "schema_version": 2, "symbols": ["O", "H", "H"], "geometry": [0.0, 0.0, 0.0, 1.43042881, 0.0, 1.10715704, -1.43042881, 0.0, 1.10715704], "molecular_charge": 0, "molecular_multiplicity": 1, "connectivity": [[0, 1, 1], [0, 2, 1]], "fragments": [[0, 1, 2]], "name": "water / TIP3P", "atom_labels": ["O", "H1", "H2"]}'


    data = H2O.to_qcschema_json()

    if correct != data:
        print("----")
        print(data)
        print("----")

    assert data == correct


def test_acetic_acid(AceticAcid):
    """Test creating QCSchema for actic acid."""
    correct = '{"schema_name": "qcschema_molecule", "schema_version": 2, "symbols": ["C", "H", "H", "H", "C", "O", "O", "H"], "geometry": [2.0403373, 0.03420404, -0.03477096, 1.09263965, 5.92920469, 0.53157996, 1.36230356, -1.27291952, -1.48513576, 1.33263486, -0.59394092, 1.80072002, 1.07960053, 2.62653034, -0.59734243, -0.25001077, 3.23936852, -2.37500779, 1.84380578, 4.34070091, 1.11852889, 4.10524103, 0.03042459, -0.05782562], "molecular_charge": 0, "molecular_multiplicity": 1, "connectivity": [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1], [4, 5, 2], [4, 6, 1], [6, 7, 1]], "fragments": [[0, 1, 2, 3, 4, 5, 6, 7]], "name": "acetic acid / acetic acid"}'

    data = AceticAcid.to_qcschema_json()

    if correct != data:
        print("----")
        print(data)
        print("----")

    assert data == correct


def test_from_schema(configuration, AceticAcid):
    """Test creating QCSchema for actic acid."""
    correct = '{"schema_name": "qcschema_molecule", "schema_version": 2, "symbols": ["C", "H", "H", "H", "C", "O", "O", "H"], "geometry": [2.040337296754675, 0.03420404285566326, -0.03477096069304994, 1.092639645256602, 5.929204688614864, 0.5315799588562472, 1.3623035632402012, -1.2729195175455674, -1.4851357613406495, 1.3326348630836315, -0.5939409209687825, 1.80072002415257, 1.0796005349967084, 2.626530340612506, -0.5973424279931026, -0.25001076628752755, 3.2393685228275113, -2.3750077934252807, 1.843805779793958, 4.340700908257376, 1.1185288931639272, 4.10524103312944, 0.030424590606418698, -0.05782561941344175], "molecular_charge": 0, "molecular_multiplicity": 1, "connectivity": [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1], [4, 5, 2], [4, 6, 1], [6, 7, 1]], "fragments": [[0, 1, 2, 3, 4, 5, 6, 7]], "name": "acetic acid / acetic acid"}'  # noqa: E501

    configuration.from_qcschema_json(correct)

    assert configuration == AceticAcid

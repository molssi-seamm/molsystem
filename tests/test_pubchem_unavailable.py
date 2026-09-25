#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A refused PubChem request must raise PubChemUnavailableError, not 'not found'."""

import pytest

import molsystem
from molsystem import pubchem, inchi


class _Response:
    def __init__(self, status_code, text="", data=None):
        self.status_code = status_code
        self.text = text
        self._data = data or {}

    def json(self):
        return self._data


@pytest.mark.parametrize("code", [403, 429, 500, 503])
def test_from_cid_unavailable(monkeypatch, code):
    monkeypatch.setattr(pubchem.requests, "get", lambda *a, **k: _Response(code))
    db = molsystem.SystemDB(filename="file:pcu?mode=memory&cache=shared")
    conf = db.create_system("s").create_configuration("c")
    with pytest.raises(molsystem.PubChemUnavailableError, match=str(code)):
        conf.PC_from_cid(5287573)
    db.close()


def test_from_cid_not_found_is_still_runtime_error(monkeypatch):
    monkeypatch.setattr(pubchem.requests, "get", lambda *a, **k: _Response(404))
    db = molsystem.SystemDB(filename="file:pcu2?mode=memory&cache=shared")
    conf = db.create_system("s").create_configuration("c")
    with pytest.raises(RuntimeError) as e:
        conf.PC_from_cid(5287573)
    assert not isinstance(e.value, molsystem.PubChemUnavailableError)
    db.close()


def test_inchikey_unavailable(monkeypatch):
    monkeypatch.setattr(inchi.requests, "get", lambda *a, **k: _Response(503))
    db = molsystem.SystemDB(filename="file:pcu3?mode=memory&cache=shared")
    conf = db.create_system("s").create_configuration("c")
    with pytest.raises(molsystem.PubChemUnavailableError):
        conf.from_inchikey("QTBSBXVTEAMEQO-UHFFFAOYSA-N")
    db.close()

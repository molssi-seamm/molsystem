# -*- coding: utf-8 -*-

"""Open Babel's own data (SSSR, LSSR, StereoData, "MOL Chiral Flag", "OpenBabel
Symmetry Classes") is not stored as properties when a structure is read (MOPAC
session, 2026-11-19: every SDF read through SEAMM gained them); the file's own
data items are."""


def test_only_the_files_data_become_properties(configuration):
    configuration.from_smiles("c1ccccc1O")
    text = configuration.to_sdf_text()
    # A data item of the file's own
    text = text.replace("$$$$", ">  <boiling point>\n181.7\n\n$$$$")
    other = configuration.system.create_configuration(name="read")
    other.from_sdf_text(text)
    names = set(other.properties.get().keys())
    for internal in (
        "SSSR",
        "LSSR",
        "StereoData",
        "MOL Chiral Flag",
        "OpenBabel Symmetry Classes",
    ):
        assert not any(name.startswith(internal) for name in names), internal
    assert any(name.startswith("boiling point") for name in names)

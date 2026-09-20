#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenBabel mixin of the class."""

import pprint  # noqa: F401
from openbabel import openbabel
import platform

import pytest  # noqa: F401

from molsystem import openbabel_version, SystemDB

# H-Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-OH
SMILES = (
    "[NH2][C@@H](CC(=O)O)C(=O)"
    "[NH][C@@H](CCCNC(=[NH2])N)C(=O)"
    "[NH][C@@H](C(C)C)C(=O)"
    "[NH][C@@H](Cc1ccc(cc1)O)C(=O)"
    "[NH][C@@H]([C@@H](CC)C)C(=O)"
    "[NH][C@@H](CC1=CNC=[NH]1)C(=O)"
    "[N]1[C@@H](CCC1)C(=O)"
    "[NH][C@@H](Cc1ccccc1)C(=O)O"
)
sidechains = {
    "ALA_LL": "[CH3X4]",
    "ARG_LL": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH2X3]",
    "ARG_LL_DHH12": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH1X2])[NH2X3]",
    "ARG_LL_DHH22": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3]([NH2X3])=[NH1X2]",
    "ARG_LL_RNH1": "[CH2X4][CH2X4][CH2X4][NH1X3][CH0X3](=[NH2X3])[NH2X3]",
    "ASN_LL": "[CH2X4][CH0X3](=[OH0X1])[NH2X3]",
    "ASP_LL": "[CH2X4][CH0X3](=[OH0X1])[OH1X2]",
    "ASP_LL_DHD2": "[CH2X4][CH0X3](=[OH0X1])[OH0X1]",
    "CYS_LL": "[CH2X4][SH1X2]",
    "CYS_LL_DHG": "[CH2X4][SH0X1]",
    "GLN_LL": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[NH2X3]",
    "GLU_LL": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH1X2]",
    "GLU_LL_DHE2": "[CH2X4][CH2X4][CH0X3](=[OH0X1])[OH0X1]",
    "HIS_LL": "[CH2X4][CH0X3][NH1X3]=[CH1X3][NH1X3][CH1X3]",
    "HIS_LL_DHD1": "[CH2X4][CH0X3][NH0X2]=[CH1X3][NH1X3][CH1X3]",
    "HIS_LL_DHE2": "[CH2X4][CH0X3][NH1X3]=[CH1X3][NH0X2][CH1X3]",
    "ILE_LL": "[CH1X4]([CH2X4][CH3X4])[CH3X4]",
    "LEU_LL": "[CH2X4][CH1X4]([CH3X4])[CH3X4]",
    "LYS_LL": "[CH2X4][CH2X4][CH2X4][CH2X4][NH3X4]",
    "LYS_LL_DHZ3": "[CH2X4][CH2X4][CH2X4][CH2X4][NH2X3]",
    "MET_LL": "[CH2X4][CH2X4][SH0X2][CH3X4]",
    "PHE_LL": "[CH2X4][cH0X3][cH1X3][cH1X3][cH1X3][cH1X3][cH1X3]",
    "SER_LL": "[CH2X4][OH1X2]",
    "SER_LL_DHG": "[CH2X4][OH0X1]",
    "THR_LL": "[CH1X4]([OH1X2])[CH3X4]",
    "THR_LL_DHG1": "[CH1X4]([OH0X1])[CH3X4]",
    "TRP_LL": (
        "[CH2X4][CH0X3]=[CH1X3][NH1X3][CH0X3]=[CH0X3][CH1X3]=[CH1X3]" "[CH1X3]=[CH1X3]"
    ),
    "TRP_LL_DHE1": (
        "[CH2X4][CH0X3]=[CH1X3][NH0X2][CH0X3]=[CH0X3][CH1X3]=[CH1X3]" "[CH1X3]=[CH1X3]"
    ),
    "TYR_LL": "[CH2X4][cH0X3]1[cH1X3][cH1X3][cH0X3]([OH1X2])[cH1X3][cH1X3]1",
    "TYR_LL_DHH": ("[CH2X4][CH0X3]=[CH1X3][CH1X3]=[CH0X3]([CH1X3]=[CH1X3])[OH0X1]"),
    "VAL_LL": "[CH1X4]([CH3X4])[CH3X4]",
}

full_smarts = {
    "GLY_LL": "[NH1X2][CH2X4][CH0X2]=[OH0X1]",
    "PRO_LL": "[NH0X2][CH1X4]([CH0X2]=[OH0X1])[CH2X4][CH2X4][CH2X4]",
}

known_input_formats = {}
known_input_formats["Darwin"] = (
    "abinit -- ABINIT Output Format",
    "acesout -- ACES output format",
    "acr -- ACR format",
    "adfband -- ADF Band output format",
    "adfdftb -- ADF DFTB output format",
    "adfout -- ADF output format",
    "alc -- Alchemy format",
    "aoforce -- Turbomole AOFORCE output format",
    "arc -- Accelrys/MSI Biosym/Insight II CAR format",
    "axsf -- XCrySDen Structure Format",
    "bgf -- MSI BGF format",
    "box -- Dock 3.5 Box format",
    "bs -- Ball and Stick format",
    "c09out -- Crystal 09 output format",
    "c3d1 -- Chem3D Cartesian 1 format",
    "c3d2 -- Chem3D Cartesian 2 format",
    "caccrt -- Cacao Cartesian format",
    "can -- Canonical SMILES format",
    "car -- Accelrys/MSI Biosym/Insight II CAR format",
    "castep -- CASTEP format",
    "ccc -- CCC format",
    "cdjson -- ChemDoodle JSON",
    "cdx -- ChemDraw binary format",
    "cdxml -- ChemDraw CDXML format",
    "cif -- Crystallographic Information File",
    "ck -- ChemKin format",
    "cml -- Chemical Markup Language",
    "cmlr -- CML Reaction format",
    "cof -- Culgi object file format",
    "CONFIG -- DL-POLY CONFIG",
    "CONTCAR -- VASP format",
    "CONTFF -- MDFF format",
    "crk2d -- Chemical Resource Kit diagram(2D)",
    "crk3d -- Chemical Resource Kit 3D format",
    "ct -- ChemDraw Connection Table format",
    "cub -- Gaussian cube format",
    "cube -- Gaussian cube format",
    "dallog -- DALTON output format",
    "dalmol -- DALTON input format",
    "dat -- Generic Output file format",
    "dmol -- DMol3 coordinates format",
    "dx -- OpenDX cube format for APBS",
    "ent -- Protein Data Bank format",
    "exyz -- Extended XYZ cartesian coordinates format",
    "fa -- FASTA format",
    "fasta -- FASTA format",
    "fch -- Gaussian formatted checkpoint file format",
    "fchk -- Gaussian formatted checkpoint file format",
    "fck -- Gaussian formatted checkpoint file format",
    "feat -- Feature format",
    "fhiaims -- FHIaims XYZ format",
    "fract -- Free Form Fractional format",
    "fs -- Fastsearch format",
    "fsa -- FASTA format",
    "g03 -- Gaussian Output",
    "g09 -- Gaussian Output",
    "g16 -- Gaussian Output",
    "g92 -- Gaussian Output",
    "g94 -- Gaussian Output",
    "g98 -- Gaussian Output",
    "gal -- Gaussian Output",
    "gam -- GAMESS Output",
    "gamess -- GAMESS Output",
    "gamin -- GAMESS Input",
    "gamout -- GAMESS Output",
    "got -- GULP format",
    "gpr -- Ghemical format",
    "gro -- GRO format",
    "gukin -- GAMESS-UK Input",
    "gukout -- GAMESS-UK Output",
    "gzmat -- Gaussian Z-Matrix Input",
    "hin -- HyperChem HIN format",
    "HISTORY -- DL-POLY HISTORY",
    "inchi -- InChI format",
    "inp -- GAMESS Input",
    "ins -- ShelX format",
    "jin -- Jaguar input format",
    "jout -- Jaguar output format",
    "log -- Generic Output file format",
    "lpmd -- LPMD format",
    "mcdl -- MCDL format",
    "mcif -- Macromolecular Crystallographic Info",
    "MDFF -- MDFF format",
    "mdl -- MDL MOL format",
    "ml2 -- Sybyl Mol2 format",
    "mmcif -- Macromolecular Crystallographic Info",
    "mmd -- MacroModel format",
    "mmod -- MacroModel format",
    "mol -- MDL MOL format",
    "mol2 -- Sybyl Mol2 format",
    "mold -- Molden format",
    "molden -- Molden format",
    "molf -- Molden format",
    "moo -- MOPAC Output format",
    "mop -- MOPAC Cartesian format",
    "mopcrt -- MOPAC Cartesian format",
    "mopin -- MOPAC Internal",
    "mopout -- MOPAC Output format",
    "mpc -- MOPAC Cartesian format",
    "mpo -- Molpro output format",
    "mpqc -- MPQC output format",
    "mrv -- Chemical Markup Language",
    "msi -- Accelrys/MSI Cerius II MSI format",
    "nwo -- NWChem output format",
    "orca -- ORCA output format",
    "out -- Generic Output file format",
    "outmol -- DMol3 coordinates format",
    "output -- Generic Output file format",
    "pc -- PubChem format",
    "pcjson -- PubChem JSON",
    "pcm -- PCModel Format",
    "pdb -- Protein Data Bank format",
    "pdbqt -- AutoDock PDBQT format",
    "png -- PNG 2D depiction",
    "pos -- POS cartesian coordinates format",
    "POSCAR -- VASP format",
    "POSFF -- MDFF format",
    "pqr -- PQR format",
    "pqs -- Parallel Quantum Solutions format",
    "prep -- Amber Prep format",
    "pwscf -- PWscf format",
    "qcout -- Q-Chem output format",
    "res -- ShelX format",
    "rsmi -- Reaction SMILES format",
    "rxn -- MDL RXN format",
    "sd -- MDL MOL format",
    "sdf -- MDL MOL format",
    "siesta -- SIESTA format",
    "smi -- SMILES format",
    "smiles -- SMILES format",
    "smy -- SMILES format using Smiley parser",
    "sy2 -- Sybyl Mol2 format",
    "t41 -- ADF TAPE41 format",
    "tdd -- Thermo format",
    "text -- Read and write raw text",
    "therm -- Thermo format",
    "tmol -- TurboMole Coordinate format",
    "txt -- Title format",
    "txyz -- Tinker XYZ format",
    "unixyz -- UniChem XYZ format",
    "VASP -- VASP format",
    "vmol -- ViewMol format",
    "wln -- Wiswesser Line Notation",
    "xml -- General XML format",
    "xsf -- XCrySDen Structure Format",
    "xyz -- XYZ cartesian coordinates format",
    "yob -- YASARA.org YOB format",
)
known_input_formats["Linux"] = (
    "abinit -- ABINIT Output Format",
    "acesout -- ACES output format",
    "acr -- ACR format",
    "adfband -- ADF Band output format",
    "adfdftb -- ADF DFTB output format",
    "adfout -- ADF output format",
    "alc -- Alchemy format",
    "aoforce -- Turbomole AOFORCE output format",
    "arc -- Accelrys/MSI Biosym/Insight II CAR format",
    "axsf -- XCrySDen Structure Format",
    "bgf -- MSI BGF format",
    "box -- Dock 3.5 Box format",
    "bs -- Ball and Stick format",
    "c09out -- Crystal 09 output format",
    "c3d1 -- Chem3D Cartesian 1 format",
    "c3d2 -- Chem3D Cartesian 2 format",
    "caccrt -- Cacao Cartesian format",
    "can -- Canonical SMILES format",
    "car -- Accelrys/MSI Biosym/Insight II CAR format",
    "castep -- CASTEP format",
    "ccc -- CCC format",
    "cdjson -- ChemDoodle JSON",
    "cdx -- ChemDraw binary format",
    "cdxml -- ChemDraw CDXML format",
    "cif -- Crystallographic Information File",
    "ck -- ChemKin format",
    "cml -- Chemical Markup Language",
    "cmlr -- CML Reaction format",
    "cof -- Culgi object file format",
    "CONFIG -- DL-POLY CONFIG",
    "CONTCAR -- VASP format",
    "CONTFF -- MDFF format",
    "crk2d -- Chemical Resource Kit diagram(2D)",
    "crk3d -- Chemical Resource Kit 3D format",
    "ct -- ChemDraw Connection Table format",
    "cub -- Gaussian cube format",
    "cube -- Gaussian cube format",
    "dallog -- DALTON output format",
    "dalmol -- DALTON input format",
    "dat -- Generic Output file format",
    "dmol -- DMol3 coordinates format",
    "dx -- OpenDX cube format for APBS",
    "ent -- Protein Data Bank format",
    "exyz -- Extended XYZ cartesian coordinates format",
    "fa -- FASTA format",
    "fasta -- FASTA format",
    "fch -- Gaussian formatted checkpoint file format",
    "fchk -- Gaussian formatted checkpoint file format",
    "fck -- Gaussian formatted checkpoint file format",
    "feat -- Feature format",
    "fhiaims -- FHIaims XYZ format",
    "fract -- Free Form Fractional format",
    "fs -- Fastsearch format",
    "fsa -- FASTA format",
    "g03 -- Gaussian Output",
    "g09 -- Gaussian Output",
    "g16 -- Gaussian Output",
    "g92 -- Gaussian Output",
    "g94 -- Gaussian Output",
    "g98 -- Gaussian Output",
    "gal -- Gaussian Output",
    "gam -- GAMESS Output",
    "gamess -- GAMESS Output",
    "gamin -- GAMESS Input",
    "gamout -- GAMESS Output",
    "got -- GULP format",
    "gpr -- Ghemical format",
    "gro -- GRO format",
    "gukin -- GAMESS-UK Input",
    "gukout -- GAMESS-UK Output",
    "gzmat -- Gaussian Z-Matrix Input",
    "hin -- HyperChem HIN format",
    "HISTORY -- DL-POLY HISTORY",
    "inchi -- InChI format",
    "inp -- GAMESS Input",
    "ins -- ShelX format",
    "jin -- Jaguar input format",
    "jout -- Jaguar output format",
    "log -- Generic Output file format",
    "lpmd -- LPMD format",
    "mcdl -- MCDL format",
    "mcif -- Macromolecular Crystallographic Info",
    "MDFF -- MDFF format",
    "mdl -- MDL MOL format",
    "ml2 -- Sybyl Mol2 format",
    "mmcif -- Macromolecular Crystallographic Info",
    "mmd -- MacroModel format",
    "mmod -- MacroModel format",
    "mol -- MDL MOL format",
    "mol2 -- Sybyl Mol2 format",
    "mold -- Molden format",
    "molden -- Molden format",
    "molf -- Molden format",
    "moo -- MOPAC Output format",
    "mop -- MOPAC Cartesian format",
    "mopcrt -- MOPAC Cartesian format",
    "mopin -- MOPAC Internal",
    "mopout -- MOPAC Output format",
    "mpc -- MOPAC Cartesian format",
    "mpo -- Molpro output format",
    "mpqc -- MPQC output format",
    "mrv -- Chemical Markup Language",
    "msi -- Accelrys/MSI Cerius II MSI format",
    "nwo -- NWChem output format",
    "orca -- ORCA output format",
    "out -- Generic Output file format",
    "outmol -- DMol3 coordinates format",
    "output -- Generic Output file format",
    "pc -- PubChem format",
    "pcjson -- PubChem JSON",
    "pcm -- PCModel Format",
    "pdb -- Protein Data Bank format",
    "pdbqt -- AutoDock PDBQT format",
    "png -- PNG 2D depiction",
    "pos -- POS cartesian coordinates format",
    "POSCAR -- VASP format",
    "POSFF -- MDFF format",
    "pqr -- PQR format",
    "pqs -- Parallel Quantum Solutions format",
    "prep -- Amber Prep format",
    "pwscf -- PWscf format",
    "qcout -- Q-Chem output format",
    "res -- ShelX format",
    "rsmi -- Reaction SMILES format",
    "rxn -- MDL RXN format",
    "sd -- MDL MOL format",
    "sdf -- MDL MOL format",
    "siesta -- SIESTA format",
    "smi -- SMILES format",
    "smiles -- SMILES format",
    "smy -- SMILES format using Smiley parser",
    "sy2 -- Sybyl Mol2 format",
    "t41 -- ADF TAPE41 format",
    "tdd -- Thermo format",
    "text -- Read and write raw text",
    "therm -- Thermo format",
    "tmol -- TurboMole Coordinate format",
    "txt -- Title format",
    "txyz -- Tinker XYZ format",
    "unixyz -- UniChem XYZ format",
    "VASP -- VASP format",
    "vmol -- ViewMol format",
    "wln -- Wiswesser Line Notation",
    "xml -- General XML format",
    "xsf -- XCrySDen Structure Format",
    "xtc -- XTC format",
    "xyz -- XYZ cartesian coordinates format",
    "yob -- YASARA.org YOB format",
)
known_input_formats["Windows"] = known_input_formats["Darwin"]

copper_sdf = """SEAMM=default/FCC Copper
 OpenBabel03252515063D

  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    1.8075    1.8075    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    1.8075    0.0000    1.8075 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.8075    1.8075 Cu  0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <SEAMM|net charge|int|>
0

>  <SEAMM|spin multiplicity|int|>
1

>  <SEAMM|XYZ|json|>
[
    [0, 0, 0],
    [1.807455, 1.807455, 0],
    [1.807455, 0, 1.807455],
    [0, 1.807455, 1.807455]
]

>  <SEAMM|cell|json|>
[3.61491, 3.61491, 3.61491, 90, 90, 90]

>  <SEAMM|system name|str|>
default

>  <SEAMM|configuration name|str|>
FCC Copper

$$$$"""

copper_sdf_2 = """SEAMM=default/FCC Copper
 OpenBabel03252519453D

  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    1.8075    1.8075    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    1.8075    3.6149    1.8075 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    3.6149    1.8075    1.8075 Cu  0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <SEAMM|net charge|int|>
0

>  <SEAMM|spin multiplicity|int|>
1

>  <SEAMM|XYZ|json|>
[
    [0, 0, 0],
    [1.807455, 1.807455, 0],
    [1.807455, 0, 1.807455],
    [0, 1.807455, 1.807455]
]

>  <SEAMM|cell|json|>
[3.61491, 3.61491, 3.61491, 90, 90, 90]

>  <SEAMM|system name|str|>
default

>  <SEAMM|configuration name|str|>
FCC Copper

$$$$"""


def test_version():
    """Test the version number of Open Babel."""
    openbabel_version()


def test_substructure(CH3COOH_3H2O):
    """Test the finding substructures in the configuration."""
    answer1 = [(5, 6, 7)]
    answer2 = [(5, 6, 7, 8)]
    answer3 = [(9,), (12,), (15,)]
    answer4 = [(10, 9, 11), (13, 12, 14), (16, 15, 17)]
    answer5 = [(9, 10, 11), (12, 13, 14), (15, 16, 17)]

    configuration = CH3COOH_3H2O

    # Just the C and O of the carboxyl, not to H on O
    result = configuration.find_substructures("C(=O)O")
    if result != answer1:
        pprint.pprint(result)
    assert result == answer1

    # All four atoms of the carboxyl group
    result = configuration.find_substructures("C(=O)[O][H]")
    if result != answer2:
        pprint.pprint(result)
    assert result == answer2

    # The Oxygen atoms of the waters
    result = configuration.find_substructures("[OH2]")
    if result != answer3:
        pprint.pprint(result)
    assert result == answer3

    # All the atoms in waters, order H-O-H
    result = configuration.find_substructures("[H][O][H]")
    if result != answer4:
        pprint.pprint(result)
    assert result == answer4

    # All the atoms in waters, order O-H-H
    result = configuration.find_substructures("[O]([H])[H]")
    if result != answer5:
        pprint.pprint(result)
    assert result == answer5


def test_substructure_ordering(disordered):
    """Test the ordering of atoms in subsets."""
    answer1 = "CC(=O)O"
    answer2 = [(1, 5, 6, 7), (16, 12, 11, 10)]
    answer3 = "C([H])([H])([H])C(=O)O[H]"
    answer4 = [(1, 2, 3, 4, 5, 6, 7, 8), (16, 15, 14, 13, 12, 11, 10, 9)]

    configuration = disordered
    templates = configuration.create_molecule_templates(create_subsets=False)

    # Without hydrogens
    smiles = templates[0].smiles
    if smiles != answer1:
        pprint.pprint(smiles)
    assert smiles == answer1

    result = configuration.find_substructures(smiles)
    if result != answer2:
        pprint.pprint(result)
    assert result == answer2

    # With hydrogens
    smiles = templates[0].to_smiles(hydrogens=True, flavor="openbabel")
    if smiles != answer3:
        pprint.pprint(smiles)
    assert smiles == answer3

    result = configuration.find_substructures(smiles)
    if result != answer4:
        pprint.pprint(result)
    assert result == answer4


def test_all_residue_search(configuration):
    """Testing locating residues in a peptide."""
    residues = {
        "ARG_LL": [(9, 10, 11, 12, 13, 14, 15, 17, 16, 18, 19)],
        "ARG_LL_RNH1": [(9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)],
        "HIS_LL": [(47, 48, 49, 50, 54, 53, 52, 51, 55, 56)],
        "ILE_LL": [(39, 40, 41, 42, 43, 44, 45, 46)],
        "PHE_LL": [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74)],
        "TYR_LL": [(27, 28, 29, 30, 31, 32, 33, 36, 34, 35, 37, 38)],
        "VAL_LL": [(20, 21, 22, 23, 24, 25, 26)],
    }
    n_terminal = {
        "ASP_LL": [(1, 2, 3, 4, 5, 6, 7, 8)],
    }
    c_terminal = {
        "PHE_LL": [(64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75)],
    }

    configuration.from_smiles(SMILES, reorient=False, flavor="openbabel")

    for name, sc in sidechains.items():
        smarts = f"[NH1X3][C@@H]({sc})[CX3]=[OX1]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in residues
        else:
            if result != residues[name]:
                print(f"'{name}': {result},")
            assert result == residues[name]

    for name, sc in sidechains.items():
        smarts = f"[NH2][C@@H]({sc})[CX3]=[OX1]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in n_terminal
        else:
            if result != n_terminal[name]:
                print(f"'{name}': {result},")
            assert result == n_terminal[name]

    for name, sc in sidechains.items():
        smarts = f"[NX3][C@@H]({sc})[CX3](=[OX1])[OH1X2]"
        result = configuration.find_substructures(smarts)
        if len(result) == 0:
            assert name not in c_terminal
        else:
            if result != c_terminal[name]:
                print(f"'{name}': {result},")
            assert result == c_terminal[name]


def test_to_OBMol(Acetate):
    """Test creating an OBMol object from a structure."""
    correct = {
        "SEAMM|XYZ|json|": (
            "[\n"
            "    [1.0797, 0.0181, -0.0184],\n"
            "    [0.5782, 3.1376, 0.2813],\n"
            "    [0.7209, -0.6736, -0.7859],\n"
            "    [0.7052, -0.3143, 0.9529],\n"
            "    [0.5713, 1.3899, -0.3161],\n"
            "    [-0.1323, 1.7142, -1.2568],\n"
            "    [0.9757, 2.297, 0.5919]\n"
            "]"
        ),
        "SEAMM|configuration name|str|": "acetate",
        "SEAMM|float property|float|kcal/mol": 3.14,
        "SEAMM|int property|int|": 2,
        "SEAMM|net charge|int|": -1,
        "SEAMM|spin multiplicity|int|": 1,
        "SEAMM|str property|str|": "Hi!",
        "SEAMM|system name|str|": "acetate",
    }

    mol = Acetate.to_OBMol(properties="*")

    bondorder_list = []
    for bond in openbabel.OBMolBondIter(mol):
        bondorder_list.append(bond.GetBondOrder())

    atno_list = []
    for atom in openbabel.OBMolAtomIter(mol):
        atno_list.append(atom.GetAtomicNum())

    assert Acetate.atoms.atomic_numbers == atno_list
    assert Acetate.bonds.bondorders == bondorder_list

    data = {}
    for item in mol.GetData():
        value = item.GetValue()
        try:
            value = int(value)
        except Exception:
            try:
                value = float(value)
            except Exception:
                pass
        data[item.GetAttribute()] = value
    if data != correct:
        pprint.pprint(data)
    assert data == correct


def test_from_OBMol(configuration):
    """Test creating a structure from an OBMol object."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "mdl")
    obConversion.AddOption("3")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, "C=CO")

    # Add hydrogens
    mol.AddHydrogens()

    # Get coordinates for a 3-D structure
    builder = openbabel.OBBuilder()
    builder.Build(mol)

    configuration.from_OBMol(mol)

    assert configuration.n_atoms == 7
    assert configuration.bonds.bondorders == [2, 1, 1, 1, 1, 1]


# Formats whose presence depends on how Open Babel was built (libxml2 for the
# XML family; newer releases add JSON formats), so they may come and go.
build_dependent_formats = {
    "cdxml -- ChemDraw CDXML format",
    "cjson -- Chemical JSON",
    "cml -- Chemical Markup Language",
    "cmlr -- CML Reaction format",
    "ket -- Ketcher KET JSON format",
    "mrv -- Chemical Markup Language",
    "pc -- PubChem format",
    "xml -- General XML format",
}


def test_input_formats():
    """Check the list of input formats Open Babel handles.

    The list is compared to the known one for the platform, ignoring formats
    that depend on the Open Babel build.
    """
    obConversion = openbabel.OBConversion()
    formats = obConversion.GetSupportedInputFormat()
    system = platform.system()
    unexpected = set(formats) ^ set(known_input_formats[system])
    unexpected -= build_dependent_formats
    if unexpected:
        import pprint

        print(system)
        pprint.pprint(formats)
    assert unexpected == set()


def test_copper_to_sdf(copper):
    """Write a manually created configuration to molfile"""
    saved = copper.to_sdf_text()
    tmp = saved.splitlines()
    del tmp[1]
    text = "\n".join(tmp)

    correct = copper_sdf.splitlines()
    del correct[1]
    correct = "\n".join(correct)

    if text != correct:
        print(saved)
    assert text == correct


def test_copper_from_sdf(copper):
    """Read an sdf file for periodic system"""

    system = copper.system.system_db.create_system()
    configuration = system.create_configuration()
    sysname, confname = configuration.from_sdf_text(copper_sdf, properties="")
    configuration.name = confname
    configuration.system.name = sysname

    result = configuration.system.diff(copper.system)
    if result != {}:
        pprint.pprint(result)
    assert result == {}

    saved = configuration.to_sdf_text()
    tmp = saved.splitlines()
    del tmp[1]
    text = "\n".join(tmp)

    correct = copper_sdf_2.splitlines()
    del correct[1]
    correct = "\n".join(correct)

    if text != correct:
        print(saved)
    assert text == correct


def test_sdf_roundtrip_keeps_dimensionless_units(AceticAcid):
    """A property with no units must not come back with units of None.

    Writing "" and reading it back as None (NULL) made the units of the property
    disagree with its definition, which then broke the unit conversion when a
    step stored a new value for it.
    """
    configuration = AceticAcid
    configuration.properties.add(
        "statistical inefficiency#LAMMPS#oplsaa+",
        _type="float",
        units="",
        description="The statistical inefficiency.",
    )
    configuration.properties.put("statistical inefficiency#LAMMPS#oplsaa+", 3.5)
    text = configuration.to_sdf_text()

    db = SystemDB(filename="file:sdf_units_db?mode=memory&cache=shared")
    try:
        new = db.create_system(name="new").create_configuration(name="new")
        new.from_sdf_text(text)

        assert new.properties.units("statistical inefficiency#LAMMPS#oplsaa+") == ""
        assert new.properties.get("statistical inefficiency#LAMMPS#oplsaa+")[
            "statistical inefficiency#LAMMPS#oplsaa+"
        ]["value"] == pytest.approx(3.5)
    finally:
        db.close()


ACETATE_SDF = """\

 test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9834   -0.0510    0.0844 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5034   -0.0510    0.0844 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1388   -0.7604    0.9259 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2134    0.7415   -0.8559 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  CHG  1   4  -1
M  END
$$$$
"""


def test_sdf_closed_shell_has_no_radical_flag(AceticAcid):
    """A closed-shell molecule must not pick up a radical flag.

    The molecular spin multiplicity used to be put on the first atom, where
    Open Babel reads 1 as a singlet carbene rather than a closed shell, so every
    SDF came out with a spurious "RAD=1"/"M  RAD" on atom 1.
    """
    configuration = AceticAcid
    assert configuration.spin_multiplicity == 1

    text = configuration.to_sdf_text()

    assert "RAD" not in text


def test_sdf_formal_charge_stays_on_its_own_atom(configuration):
    """The molecular charge must not be written as a formal charge on atom 1."""
    configuration.from_sdf_text(ACETATE_SDF)

    assert configuration.charge == -1
    # The charge is on the 4th atom, one of the two oxygens.
    assert configuration.atoms.get_column_data("formal_charge") == [0, 0, 0, -1]

    text = configuration.to_sdf_text()

    assert "RAD" not in text
    charges = [line.strip() for line in text.splitlines() if "CHG" in line]
    assert charges == ["M  CHG  1   4  -1"]


def _water_box(configuration, n_mol=400, charged_molecule=None):
    """Fill a configuration with enough water to force the V3000 SDF format.

    Open Babel switches to V3000 above 999 atoms or bonds, so this is the format
    the production 500-molecule cells are written in.
    """
    Xs, Ys, Zs, atnos, qs = [], [], [], [], []
    geometry = ((8, (0.0, 0.0, 0.0)), (1, (0.96, 0.0, 0.0)), (1, (-0.24, 0.93, 0.0)))
    for molecule in range(n_mol):
        x0 = molecule * 5.0
        for atno, (dx, dy, dz) in geometry:
            atnos.append(atno)
            Xs.append(x0 + dx)
            Ys.append(dy)
            Zs.append(dz)
            qs.append(0)

    if charged_molecule is not None:
        # Make one molecule a hydroxide: the charge sits on its oxygen.
        qs[3 * charged_molecule] = -1
        configuration.atoms.add_attribute("formal_charge", coltype="int", default=0)
        ids = configuration.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos, formal_charge=qs)
        configuration.charge = -1
    else:
        ids = configuration.atoms.append(x=Xs, y=Ys, z=Zs, atno=atnos)

    Is = [ids[3 * m] for m in range(n_mol)] * 2
    Js = [ids[3 * m + 1] for m in range(n_mol)]
    Js += [ids[3 * m + 2] for m in range(n_mol)]
    configuration.bonds.append(i=Is, j=Js, bondorder=[1] * (2 * n_mol))


def test_v3000_sdf_closed_shell_has_no_radical_flag(configuration):
    """The V3000 writer showed the same atom-1 radical flag as V2000.

    This is the form it was found in: "M  V30 1 O ... RAD=1" in the 500-molecule
    liquid cells, which are large enough that Open Babel picks V3000.
    """
    _water_box(configuration)

    text = configuration.to_sdf_text()

    assert "V3000" in text
    assert "RAD" not in text


def test_v3000_sdf_formal_charge_stays_on_its_own_atom(configuration):
    """In V3000 too, the molecular charge must not land on atom 1."""
    _water_box(configuration, charged_molecule=99)

    assert configuration.charge == -1

    text = configuration.to_sdf_text()

    assert "V3000" in text
    assert "RAD" not in text
    charges = [line.strip() for line in text.splitlines() if "CHG=" in line]
    # The oxygen of the 100th water, i.e. the 298th atom -- not the first.
    assert charges == ["M  V30 298 O 495 0 0 0 CHG=-1"]

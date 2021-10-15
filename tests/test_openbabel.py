#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenBabel mixin of the class."""

import pprint  # noqa: F401
from openbabel import openbabel
import platform

import pytest  # noqa: F401

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
    smiles = templates[0].to_smiles(hydrogens=True)
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

    configuration.from_smiles(SMILES)

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


def test_to_OBMol(configuration):
    """Test creating an OBMol object from a structure."""
    mol = configuration.to_OBMol()

    bondorder_list = []
    for bond in openbabel.OBMolBondIter(mol):
        bondorder_list.append(bond.GetBondOrder())

    atno_list = []
    for atno in openbabel.OBMolAtomIter(mol):
        atno_list.append(mol.GetAtmoicNum())

    assert configuration.atoms.atomic_numbers == atno_list
    assert configuration.bonds.bondorders == bondorder_list
    # assert configuration.n_atoms == mol.NumAtoms()
    # assert configuration.n_bonds == mol.NumBonds()


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


def test_input_formats():
    """Check the list of input formats Open Babel handles"""
    obConversion = openbabel.OBConversion()
    formats = obConversion.GetSupportedInputFormat()
    system = platform.system()
    if formats != known_input_formats[system]:
        import pprint

        print(system)
        pprint.pprint(formats)
    assert formats == known_input_formats[system]

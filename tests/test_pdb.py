#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

"""Tests for handling PDB files."""

ruler = """\
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
#-----+++++ ----#--- #----#   --------########--------######------          --##
"""  # noqa: E501

# https://files.rcsb.org/ligands/download/HEM_ideal.pdb

heme = """\
ATOM      1  CHA HEM A   1      -2.161  -0.125   0.490  1.00 10.00           C
ATOM      2  CHB HEM A   1       1.458  -3.419   0.306  1.00 10.00           C
ATOM      3  CHC HEM A   1       4.701   0.169  -0.069  1.00 10.00           C
ATOM      4  CHD HEM A   1       1.075   3.460   0.018  1.00 10.00           C
ATOM      5  C1A HEM A   1      -1.436  -1.305   0.380  1.00 10.00           C
ATOM      6  C2A HEM A   1      -2.015  -2.587   0.320  1.00 10.00           C
ATOM      7  C3A HEM A   1      -1.009  -3.500   0.270  1.00 10.00           C
ATOM      8  C4A HEM A   1       0.216  -2.803   0.298  1.00 10.00           C
ATOM      9  CMA HEM A   1      -1.175  -4.996   0.197  1.00 10.00           C
ATOM     10  CAA HEM A   1      -3.490  -2.893   0.314  1.00 10.00           C
ATOM     11  CBA HEM A   1      -3.998  -2.926  -1.129  1.00 10.00           C
ATOM     12  CGA HEM A   1      -5.473  -3.232  -1.136  1.00 10.00           C
ATOM     13  O1A HEM A   1      -6.059  -3.405  -0.094  1.00 10.00           O
ATOM     14  O2A HEM A   1      -6.137  -3.311  -2.300  1.00 10.00           O
ATOM     15  C1B HEM A   1       2.664  -2.707   0.308  1.00 10.00           C
ATOM     16  C2B HEM A   1       3.937  -3.328   0.418  1.00 10.00           C
ATOM     17  C3B HEM A   1       4.874  -2.341   0.314  1.00 10.00           C
ATOM     18  C4B HEM A   1       4.117  -1.079   0.139  1.00 10.00           C
ATOM     19  CMB HEM A   1       4.203  -4.798   0.613  1.00 10.00           C
ATOM     20  CAB HEM A   1       6.339  -2.497   0.365  1.00 10.00           C
ATOM     21  CBB HEM A   1       6.935  -3.419  -0.385  1.00 10.00           C
ATOM     22  C1C HEM A   1       3.964   1.345  -0.174  1.00 10.00           C
ATOM     23  C2C HEM A   1       4.531   2.601  -0.445  1.00 10.00           C
ATOM     24  C3C HEM A   1       3.510   3.536  -0.437  1.00 10.00           C
ATOM     25  C4C HEM A   1       2.304   2.846  -0.139  1.00 10.00           C
ATOM     26  CMC HEM A   1       5.991   2.880  -0.697  1.00 10.00           C
ATOM     27  CAC HEM A   1       3.649   4.981  -0.692  1.00 10.00           C
ATOM     28  CBC HEM A   1       4.201   5.407  -1.823  1.00 10.00           C
ATOM     29  C1D HEM A   1      -0.102   2.753   0.298  1.00 10.00           C
ATOM     30  C2D HEM A   1      -1.382   3.388   0.641  1.00 10.00           C
ATOM     31  C3D HEM A   1      -2.283   2.389   0.774  1.00 10.00           C
ATOM     32  C4D HEM A   1      -1.561   1.137   0.511  1.00 10.00           C
ATOM     33  CMD HEM A   1      -1.639   4.863   0.811  1.00 10.00           C
ATOM     34  CAD HEM A   1      -3.741   2.532   1.123  1.00 10.00           C
ATOM     35  CBD HEM A   1      -4.573   2.563  -0.160  1.00 10.00           C
ATOM     36  CGD HEM A   1      -6.032   2.706   0.189  1.00 10.00           C
ATOM     37  O1D HEM A   1      -6.372   2.776   1.347  1.00 10.00           O
ATOM     38  O2D HEM A   1      -6.954   2.755  -0.785  1.00 10.00           O
ATOM     39  NA  HEM A   1      -0.068  -1.456   0.321  1.00 10.00           N
ATOM     40  NB  HEM A   1       2.820  -1.386   0.207  1.00 10.00           N
ATOM     41  NC  HEM A   1       2.604   1.506  -0.033  1.00 10.00           N
ATOM     42  ND  HEM A   1      -0.276   1.431   0.298  1.00 10.00           N
ATOM     43 FE   HEM A   1       1.010   0.157  -0.060  1.00 10.00          FE
ATOM     44  HHB HEM A   1       1.498  -4.508   0.309  1.00 10.00           H
ATOM     45  HHC HEM A   1       5.786   0.229  -0.153  1.00 10.00           H
ATOM     46  HHD HEM A   1       1.018   4.543  -0.083  1.00 10.00           H
ATOM     47  HMA HEM A   1      -1.220  -5.306  -0.847  1.00 10.00           H
ATOM     48 HMAA HEM A   1      -0.328  -5.480   0.683  1.00 10.00           H
ATOM     49 HMAB HEM A   1      -2.097  -5.285   0.702  1.00 10.00           H
ATOM     50  HAA HEM A   1      -3.662  -3.862   0.782  1.00 10.00           H
ATOM     51 HAAA HEM A   1      -4.024  -2.121   0.869  1.00 10.00           H
ATOM     52  HBA HEM A   1      -3.825  -1.956  -1.597  1.00 10.00           H
ATOM     53 HBAA HEM A   1      -3.464  -3.697  -1.684  1.00 10.00           H
ATOM     54  HMB HEM A   1       3.256  -5.336   0.660  1.00 10.00           H
ATOM     55 HMBA HEM A   1       4.794  -5.175  -0.222  1.00 10.00           H
ATOM     56 HMBB HEM A   1       4.752  -4.948   1.543  1.00 10.00           H
ATOM     57  HAB HEM A   1       6.927  -1.863   1.011  1.00 10.00           H
ATOM     58  HBB HEM A   1       7.994  -3.600  -0.277  1.00 10.00           H
ATOM     59 HBBA HEM A   1       6.360  -3.987  -1.102  1.00 10.00           H
ATOM     60  HMC HEM A   1       6.554   1.949  -0.639  1.00 10.00           H
ATOM     61 HMCA HEM A   1       6.110   3.316  -1.689  1.00 10.00           H
ATOM     62 HMCB HEM A   1       6.362   3.578   0.053  1.00 10.00           H
ATOM     63  HAC HEM A   1       3.303   5.694   0.042  1.00 10.00           H
ATOM     64  HBC HEM A   1       4.614   4.696  -2.523  1.00 10.00           H
ATOM     65 HBCA HEM A   1       4.235   6.464  -2.043  1.00 10.00           H
ATOM     66  HMD HEM A   1      -0.715   5.415   0.639  1.00 10.00           H
ATOM     67 HMDA HEM A   1      -2.394   5.185   0.094  1.00 10.00           H
ATOM     68 HMDB HEM A   1      -1.994   5.055   1.824  1.00 10.00           H
ATOM     69  HAD HEM A   1      -4.052   1.687   1.738  1.00 10.00           H
ATOM     70 HADA HEM A   1      -3.893   3.459   1.677  1.00 10.00           H
ATOM     71  HBD HEM A   1      -4.262   3.408  -0.775  1.00 10.00           H
ATOM     72 HBDA HEM A   1      -4.421   1.636  -0.714  1.00 10.00           H
ATOM     73  H2A HEM A   1      -7.082  -3.510  -2.254  1.00 10.00           H
ATOM     74  H2D HEM A   1      -7.877   2.847  -0.512  1.00 10.00           H
ATOM     75  HHA HEM A   1      -3.246  -0.188   0.567  1.00 10.00           H
CONECT    1    5   32   75
CONECT    2    8   15   44
CONECT    3   18   22   45
CONECT    4   25   29   46
CONECT    5    1    6   39
CONECT    6    5    7   10
CONECT    7    6    8    9
CONECT    8    2    7   39
CONECT    9    7   47   48   49
CONECT   10    6   11   50   51
CONECT   11   10   12   52   53
CONECT   12   11   13   14
CONECT   13   12
CONECT   14   12   73
CONECT   15    2   16   40
CONECT   16   15   17   19
CONECT   17   16   18   20
CONECT   18    3   17   40
CONECT   19   16   54   55   56
CONECT   20   17   21   57
CONECT   21   20   58   59
CONECT   22    3   23   41
CONECT   23   22   24   26
CONECT   24   23   25   27
CONECT   25    4   24   41
CONECT   26   23   60   61   62
CONECT   27   24   28   63
CONECT   28   27   64   65
CONECT   29    4   30   42
CONECT   30   29   31   33
CONECT   31   30   32   34
CONECT   32    1   31   42
CONECT   33   30   66   67   68
CONECT   34   31   35   69   70
CONECT   35   34   36   71   72
CONECT   36   35   37   38
CONECT   37   36
CONECT   38   36   74
CONECT   39    5    8   43
CONECT   40   15   18   43
CONECT   41   22   25   43
CONECT   42   29   32   43
CONECT   43   39   40   41   42
CONECT   44    2
CONECT   45    3
CONECT   46    4
CONECT   47    9
CONECT   48    9
CONECT   49    9
CONECT   50   10
CONECT   51   10
CONECT   52   11
CONECT   53   11
CONECT   54   19
CONECT   55   19
CONECT   56   19
CONECT   57   20
CONECT   58   21
CONECT   59   21
CONECT   60   26
CONECT   61   26
CONECT   62   26
CONECT   63   27
CONECT   64   28
CONECT   65   28
CONECT   66   33
CONECT   67   33
CONECT   68   33
CONECT   69   34
CONECT   70   34
CONECT   71   35
CONECT   72   35
CONECT   73   14
CONECT   74   38
CONECT   75    1
END"""


def test_to_pdb_text(AceticAcid):
    """Write a manually created system to pdb"""
    correct = """\
COMPND    UNNAMED
AUTHOR    MolSSI SEAMM at 0915201350
ATOM      1  C   UNK A   1       1.080   0.018  -0.018  1.00  0.00           C  
ATOM      2  H   UNK A   1       0.578   3.138   0.281  1.00  0.00           H  
ATOM      3  H   UNK A   1       0.721  -0.674  -0.786  1.00  0.00           H  
ATOM      4  H   UNK A   1       0.705  -0.314   0.953  1.00  0.00           H  
ATOM      5  C   UNK A   1       0.571   1.390  -0.316  1.00  0.00           C  
ATOM      6  O   UNK A   1      -0.132   1.714  -1.257  1.00  0.00           O  
ATOM      7  O   UNK A   1       0.976   2.297   0.592  1.00  0.00           O  
ATOM      8  H   UNK A   1       2.172   0.016  -0.031  1.00  0.00           H  
CONECT    1    2    3    4    5
CONECT    2    1
CONECT    3    1
CONECT    4    1
CONECT    5    1    6    7
CONECT    6    5
CONECT    7    5    8
CONECT    8    7
MASTER        0    0    0    0    0    0    0    0    8    0    8
END"""  # noqa: E501, W291

    text = AceticAcid.to_pdb_text()

    tmp = correct.splitlines()
    del tmp[1]
    tmp_correct = "\n".join(tmp)

    tmp = text.splitlines()
    del tmp[1]
    tmp_text = "\n".join(tmp)

    if tmp_text != tmp_correct:
        print(text)

    assert tmp_text == tmp_correct


def test_from_pdb_text(configuration):
    """Test reading a standard PDB file"""
    configuration.from_pdb_text(heme)
    assert configuration.n_atoms == 75


def test_from_to(configuration):
    """Test that getting from then putting to a string works"""
    configuration.from_pdb_text(heme)
    text = configuration.to_pdb_text()

    # Strip the first two and next to last line from SEAMM's
    # representation since the PDB example ignores. Also any trailing blanks
    tmp = text.splitlines()
    del tmp[0:2]
    del tmp[-2]
    tmp2 = []
    for line in tmp:
        tmp2.append(line.rstrip())
    tmp_text = "\n".join(tmp2)

    if tmp_text != heme:
        print("diffs")
        i = 0
        for old, new in zip(heme.splitlines(), tmp2):
            i += 1
            if old != new:
                print(f"{i:3} {old}\n    {new}")

    assert tmp_text == heme

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest  # noqa: F401

"""Tests for handling Molfiles."""

text1 = """\
L-Alanine
GSMACCS-II07189510252D 1 0.00366 0.00000 0
Figure 1, J. Chem. Inf. Comput. Sci., Vol 32, No. 3., 1992
0 0 0 0 0 999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -0.6622 0.5342 0 0 CFG=2
M  V30 2 C 0.6622 -0.3 0 0
M  V30 3 C -0.7207 2.0817 0 0 MASS=13
M  V30 4 N -1.8622 -0.3695 0 0 CHG=1
M  V30 5 O 0.622 -1.8037 0 0
M  V30 6 O 1.9464 0.4244 0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3 CFG=1
M  V30 3 1 1 4
M  V30 4 2 2 5
M  V30 5 1 2 6
M  V30 END BOND
M  V30 END CTAB
M  END
"""

atoms1 = """\
   atno  formal_charge  configuration       x       y    z
1     6              0              1 -0.6622  0.5342  0.0
2     6              0              1  0.6622 -0.3000  0.0
3     6              0              1 -0.7207  2.0817  0.0
4     7              1              1 -1.8622 -0.3695  0.0
5     8              0              1  0.6220 -1.8037  0.0
6     8             -1              1  1.9464  0.4244  0.0"""

bonds1 = """\
   i  j  bondorder symop1 symop2
1  1  2          1      .      .
2  1  3          1      .      .
3  1  4          1      .      .
4  2  5          2      .      .
5  2  6          1      .      ."""

text2 = """\
L-Alanine
PSSEAMM_WF09112015563D
Exported from SEAMM
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.6622 0.5342 0.0 0
M  V30 2 C 0.6622 -0.3 0.0 0
M  V30 3 C -0.7207 2.0817 0.0 0
M  V30 4 N -1.8622 -0.3695 0.0 0
M  V30 5 O 0.622 -1.8037 0.0 0
M  V30 6 O 1.9464 0.4244 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 2 2 5
M  V30 5 1 2 6
M  V30 END BOND
M  V30 END CTAB
M  END
"""

text3 = """\
acetic acid
PSSEAMM_WF09112018063D
Exported from SEAMM
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.0797 0.0181 -0.0184 0
M  V30 2 H 0.5782 3.1376 0.2813 0
M  V30 3 H 0.7209 -0.6736 -0.7859 0
M  V30 4 H 0.7052 -0.3143 0.9529 0
M  V30 5 C 0.5713 1.3899 -0.3161 0
M  V30 6 O -0.1323 1.7142 -1.2568 0
M  V30 7 O 0.9757 2.297 0.5919 0
M  V30 8 H 2.1724 0.0161 -0.0306 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 5 2 5 6
M  V30 6 1 5 7
M  V30 7 1 7 8
M  V30 END BOND
M  V30 END CTAB
M  END
"""


def test_from_text(configuration):
    """Create a configuration from the text of a Molfile"""
    configuration.from_molfile_text(text1)
    if str(configuration.atoms) != atoms1:
        print(configuration.atoms)
    if str(configuration.bonds) != bonds1:
        print(configuration.bonds)
    assert str(configuration.atoms) == atoms1
    assert str(configuration.bonds) == bonds1


def test_to_text(configuration):
    """Create the text of a Molfile from a configuration"""
    configuration.from_molfile_text(text1)
    # Second line contains the current date/time
    text_sv = configuration.to_molfile_text()
    text = text_sv.splitlines()
    del text[1]
    text = "\n".join(text)

    tmp = text2.splitlines()
    del tmp[1]
    tmp = "\n".join(tmp)

    if text != tmp:
        print(text_sv)
    assert text == tmp


def test_to_text2(AceticAcid):
    """Write a manually created configuration to molfile"""
    # Second line contains the current date/time
    text_sv = AceticAcid.to_molfile_text()
    text = text_sv.splitlines()
    del text[1]
    text = "\n".join(text)

    tmp = text3.splitlines()
    del tmp[1]
    tmp = "\n".join(tmp)

    if text != tmp:
        print(text_sv)
    assert text == tmp

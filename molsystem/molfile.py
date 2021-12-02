# -*- coding: utf-8 -*-

"""Functions for handling MDL molfiles"""

import logging
import time

from molsystem import elements

logger = logging.getLogger(__name__)


class MolFileMixin:
    """A mixin for handling MDL Molfiles."""

    def to_molfile_text(self, title=None, comment="Exported from SEAMM"):
        """Create the text of the Molfile from the system.

        Parameters
        ----------
        title : str = None
            The title for the structure, by default the system name.
        comment : str = 'Exported from SEAMM'
            Comment line

        Returns
        -------
        text : str
            The text of the file.
        """

        lines = []

        atoms = self.atoms
        bonds = self.bonds

        n_atoms = atoms.n_atoms
        n_bonds = bonds.n_bonds

        nsgroups = 0
        n3d = 0
        is_chiral = 0  # may need to think about this later.

        if title is None:
            lines.append(self.name)
        else:
            lines.append(title)
        date_time = time.strftime("%m%d%y%H%M")

        lines.append("PS" + "SEAMM_WF" + date_time + "3D")
        lines.append(comment)
        lines.append("  0  0  0     0  0            999 V3000")

        lines.append("M  V30 BEGIN CTAB")
        lines.append(
            "M  V30 COUNTS {} {} {} {} {}".format(
                n_atoms, n_bonds, nsgroups, n3d, is_chiral
            )
        )
        lines.append("M  V30 BEGIN ATOM")
        count = 0
        if "formal charges" in atoms:
            for row in atoms.atoms():
                count += 1
                symbol = elements.to_symbols([row["atno"]])[0]
                lines.append(
                    f"M  V30 {count} {symbol} {row['x']} {row['y']} {row['z']}"
                    " 0 CHG={row['formal charge']}"
                )
        else:
            for row in atoms.atoms():
                count += 1
                symbol = elements.to_symbols([row["atno"]])[0]
                lines.append(
                    f"M  V30 {count} {symbol} {row['x']} {row['y']} {row['z']}" " 0"
                )
        lines.append("M  V30 END ATOM")
        lines.append("M  V30 BEGIN BOND")
        count = 0
        for row in bonds.bonds():
            count += 1
            lines.append(f"M  V30 {count} {row['bondorder']} " f"{row['i']} {row['j']}")
        lines.append("M  V30 END BOND")
        lines.append("M  V30 END CTAB")
        lines.append("M  END")

        return "\n".join(lines)

    def from_molfile_text(self, data):
        """Create the system from an MDL Molfile, version 3

        Parameters
        ----------
        data : str
            The complete text of the Molfile.
        """

        self.clear()
        self.periodicity = 0

        n_molecules = 0
        lines = enumerate(data.splitlines())

        # title
        lineno, title = next(lines)
        self.name = title.strip()
        # header
        next(lines)
        # comment
        next(lines)
        lineno, line = next(lines)
        if line.split()[6] != "V3000":
            raise RuntimeError(
                f"molfile:to_seamm -- the file is not version 3: '{line}'"
            )
        for lineno, line in lines:
            logger.debug(f"{lineno}: {line}")
            if "M  END" in line:
                break
            elif "M  V30 BEGIN CTAB" in line:
                n_molecules += 1
                if n_molecules > 1:
                    raise NotImplementedError("Multiple molecules?")
            elif "M V30 END CTAB" in line:
                pass
            elif "M  V30 COUNTS" in line:
                natoms, nbonds, nsgroups, n3d, is_chiral = line.split()[3:]

                natoms = int(natoms)
                nbonds = int(nbonds)
                # not used, yet.
                # nsgroups = int(nsgroups)
                # n3d = int(n3d)
                # is_chiral = bool(is_chiral)
            elif "M  V30 BEGIN ATOM" in line:
                logger.debug("In atom table")
                xs = []
                ys = []
                zs = []
                symbols = []
                formal_charges = []
                have_formal_charges = False
                for lineno, line in lines:
                    if "M  V30 END ATOM" in line:
                        logger.debug(f"Saving {len(xs)} atoms to system")
                        if have_formal_charges and "formal_charge" not in self.atoms:
                            logger.debug("   with formal charges")
                            self.atoms.add_attribute(
                                "formal_charge", coltype="int", default=0
                            )
                            atom_ids = self.atoms.append(
                                x=xs,
                                y=ys,
                                z=zs,
                                symbol=symbols,
                                formal_charge=formal_charges,
                            )
                        else:
                            atom_ids = self.atoms.append(
                                x=xs, y=ys, z=zs, symbol=symbols
                            )
                        break
                    i, symbol, x, y, z, q = line.split()[2:8]
                    xs.append(float(x))
                    ys.append(float(y))
                    zs.append(float(z))
                    symbols.append(symbol)
                    if "CHG=" in line:
                        for tmp in line.split()[8:]:
                            if "CHG=" in tmp:
                                formal_charges.append(int(tmp[4:]))
                                have_formal_charges = True
                    else:
                        formal_charges.append(0)
            elif "M  V30 BEGIN BOND" in line:
                logger.debug("In bond table")
                iatoms = []
                jatoms = []
                bondorders = []
                for lineno, line in lines:
                    if "M  V30 END BOND" in line:
                        if len(iatoms) > 0:
                            self.bonds.append(
                                i=iatoms,
                                j=jatoms,
                                bondorder=bondorders,
                            )
                        break
                    bondorder, iatom, jatom = line.split()[3:6]
                    iatoms.append(atom_ids[int(iatom) - 1])
                    jatoms.append(atom_ids[int(jatom) - 1])
                    bondorders.append(int(bondorder))

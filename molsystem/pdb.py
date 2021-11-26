# -*- coding: utf-8 -*-

"""Functions for handling PDB files

To Do
^^^^^

Need to understand more fully the PDB/mmcif format and the how to
carry the information about residues, chains, hetero groups, waters, etc.
At the moment this is ignoring much of the information, and putting residue,
chain, etc information directly on atoms.

I think we should use templates and subsets, but am not (yet) sure.

Presumably this metadata is most useful for setting up complicated
simulations.

File Format
^^^^^^^^^^^

For complete documentation, see
http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

Order of records::

    RECORD TYPE             EXISTENCE           CONDITIONS IF  OPTIONAL
    --------------------------------------------------------------------------------------
    HEADER                  Mandatory
    OBSLTE                  Optional            Mandatory in  entries that have been
                                                replaced by a newer entry.
    TITLE                   Mandatory
    SPLIT                   Optional            Mandatory when  large macromolecular
                                                complexes  are split into multiple PDB
                                                entries.
    CAVEAT                  Optional            Mandatory when there are outstanding  errors
                                                such  as chirality.
    COMPND                  Mandatory
    SOURCE                  Mandatory
    KEYWDS                  Mandatory
    EXPDTA                  Mandatory
    NUMMDL                  Optional            Mandatory for  NMR ensemble entries.
    MDLTYP                  Optional            Mandatory for  NMR minimized average
                                                Structures or when the entire  polymer
                                                chain contains C alpha or P atoms only.
    AUTHOR                  Mandatory
    REVDAT                  Mandatory
    SPRSDE                  Optional            Mandatory for a replacement entry.
    JRNL                    Optional            Mandatory for a publication describes
                                                the experiment.
    REMARK 0                Optional            Mandatory for a re-refined structure
    REMARK 1                Optional
    REMARK 2                Mandatory
    REMARK 3                Mandatory
    REMARK N                Optional            Mandatory under certain conditions.
    DBREF                   Optional            Mandatory for all polymers.
    DBREF1/DBREF2           Optional            Mandatory when certain sequence  database
                                                accession  and/or sequence numbering
                                                does  not fit preceding DBREF format.
    SEQADV                  Optional            Mandatory if sequence  conflict exists.
    SEQRES                  Mandatory           Mandatory if ATOM records exist.
    MODRES                  Optional            Mandatory if modified group exists  in the
                                                coordinates.
    HET                     Optional            Mandatory if a non-standard group other
                                                than water appears in the coordinates.
    HETNAM                  Optional            Mandatory if a non-standard group other
                                                than  water appears in the coordinates.
    HETSYN                  Optional
    FORMUL                  Optional            Mandatory if a non-standard group or
                                                water appears in the coordinates.
    HELIX                   Optional
    SHEET                   Optional
    SSBOND                  Optional            Mandatory if a  disulfide bond is present.
    LINK                    Optional            Mandatory if  non-standard residues appear
                                                in a  polymer
    CISPEP                  Optional
    SITE                    Optional
    CRYST1                  Mandatory
    ORIGX1 ORIGX2 ORIGX3    Mandatory
    SCALE1 SCALE2 SCALE3    Mandatory
    MTRIX1 MTRIX2 MTRIX3    Optional            Mandatory if  the complete asymmetric unit
                                                must  be generated from the given coordinates
                                                using non-crystallographic symmetry.
    MODEL                   Optional            Mandatory if more than one model
                                                is  present in the entry.
    ATOM                    Optional            Mandatory if standard residues exist.
    ANISOU                  Optional
    TER                     Optional            Mandatory if ATOM records exist.
    HETATM                  Optional            Mandatory if non-standard group exists.
    ENDMDL                  Optional            Mandatory if MODEL appears.
    CONECT                  Optional            Mandatory if non-standard group appears
                                                and  if LINK or SSBOND records exist.
    MASTER                  Mandatory
    END                     Mandatory

Description of HETATM records::

    COLUMNS       DATA  TYPE     FIELD         DEFINITION
    -----------------------------------------------------------------------
     1 - 6        Record name    "HETATM"
     7 - 11       Integer        serial        Atom serial number.
    13 - 16       Atom           name          Atom name.
    17            Character      altLoc        Alternate location indicator.
    18 - 20       Residue name   resName       Residue name.
    22            Character      chainID       Chain identifier.
    23 - 26       Integer        resSeq        Residue sequence number.
    27            AChar          iCode         Code for insertion of residues.
    31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
    39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
    47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
    55 - 60       Real(6.2)      occupancy     Occupancy.
    61 - 66       Real(6.2)      tempFactor    Temperature factor.
    77 - 78       LString(2)     element       Element symbol; right-justified.
    79 - 80       LString(2)     charge        Charge on the atom.

Description of CONECT records::

    COLUMNS       DATA  TYPE      FIELD        DEFINITION
    -------------------------------------------------------------------------
     1 -  6        Record name    "CONECT"
     7 - 11        Integer        serial       Atom  serial number
    12 - 16        Integer        serial       Serial number of bonded atom
    17 - 21        Integer        serial       Serial  number of bonded atom
    22 - 26        Integer        serial       Serial number of bonded atom
    27 - 31        Integer        serial       Serial number of bonded atom
"""  # noqa: E501

import collections
import logging
import time

logger = logging.getLogger(__name__)


class PDBMixin:
    """A mixin for handling PDB files."""

    def to_pdb_text(self, title=None, comment="Exported from SEAMM"):
        """Create the text of the PDB file from the system.

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

        n_atoms = atoms.n_atoms

        date_time = time.strftime("%m%d%y%H%M")

        lines.append("COMPND    UNNAMED")
        lines.append("AUTHOR    MolSSI SEAMM at " + date_time)

        # atoms

        if "resname" in atoms:
            resnames = atoms["resname"]
        else:
            resnames = ["UNK"] * n_atoms

        if "chainid" in atoms:
            chainids = atoms["chainid"]
        else:
            chainids = ["A"] * n_atoms

        if "resseq" in atoms:
            resseqs = [1 if x is None else x for x in atoms["resseq"]]
        else:
            resseqs = [1] * n_atoms

        if "occupancy" in atoms:
            occupancies = atoms["occupancy"]
        else:
            occupancies = [1.0] * n_atoms

        if "tempfactor" in atoms:
            tempfactors = atoms["tempfactor"]
        else:
            tempfactors = [0.0] * n_atoms

        if "formal_charge" in atoms:
            charges = atoms["formal_charge"]
        else:
            charges = [" "] * n_atoms

        count = 0
        symbols = atoms.symbols
        coordinates = atoms.coordinates
        if "name" in atoms:
            names = [
                symbol if name is None else name
                for name, symbol in zip(atoms["name"], symbols)
            ]
        else:
            names = symbols

        for (
            element,
            xyz,
            name,
            resname,
            chainid,
            resseq,
            occupancy,
            tempfactor,
            charge,
        ) in zip(
            symbols,
            coordinates,
            names,
            resnames,
            chainids,
            resseqs,
            occupancies,
            tempfactors,
            charges,
        ):
            count += 1
            x, y, z = xyz
            if len(element) == 1 and len(name) < 4:
                name = " " + name
            lines.append(
                f"ATOM  {count:5d} {name:4s} {resname:3s} {chainid:1s}"
                f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}"
                f"{tempfactor:6.2f}          {element.upper():>2s}"
                f"{charge:2}"
            )

        # bonds
        for i, js in enumerate(self.bonded_neighbors(as_indices=True)):
            if len(js) > 0:
                lines.append(f"CONECT{i+1:5d}" + "".join([f"{j+1:5d}" for j in js]))

        lines.append(
            "MASTER        0    0    0    0    0    0    0    0"
            f"{n_atoms:5d}    0{n_atoms:5d}"
        )
        lines.append("END")

        return "\n".join(lines)

    def from_pdb_text(self, data):
        """Create the system from a PDF file.

        Parameters
        ----------
        data : str
            The complete text of the Molfile.
        """
        # Initialize the structure

        self.clear()
        self.periodicity = 0

        if isinstance(data, list):
            lines = data
        else:
            lines = data.splitlines()

        names = []
        symbols = []
        xs = []
        ys = []
        zs = []
        resnames = []
        chainids = []
        resseqs = []
        occupancies = []
        tempfactors = []
        connections = None
        for line in lines:
            key = line[0:6].rstrip()
            if key == "HEADER":
                pass
            elif key == "OBSLTE":
                pass
            elif key == "TITLE":
                pass
            elif key == "SPLIT":
                pass
            elif key == "CAVEAT":
                pass
            elif key == "COMPND":
                pass
            elif key == "SOURCE":
                pass
            elif key == "KEYWDS":
                pass
            elif key == "EXPDTA":
                pass
            elif key == "NUMMDL":
                pass
            elif key == "MDLTYPE":
                pass
            elif key == "AUTHOR":
                pass
            elif key == "REVDAT":
                pass
            elif key == "SPRSDE":
                pass
            elif key == "JRNL":
                pass
            elif key == "REMARK":
                pass
            elif key == "DBREF":
                pass
            elif key == "DBREF1":
                pass
            elif key == "DBREF2":
                pass
            elif key == "SEQADV":
                pass
            elif key == "SEQRES":
                pass
            elif key == "MODRES":
                pass
            elif key == "HET":
                pass
            elif key == "HETNAM":
                pass
            elif key == "HETSYN":
                pass
            elif key == "FORMUL":
                pass
            elif key == "HELIX":
                pass
            elif key == "SHEET":
                pass
            elif key == "SSBOND":
                pass
            elif key == "LINK":
                pass
            elif key == "CISPEP":
                pass
            elif key == "SITE":
                pass
            elif key == "CRYST1":
                pass
            elif key == "ORIGX1":
                pass
            elif key == "ORIGX2":
                pass
            elif key == "ORIGX3":
                pass
            elif key == "SCALE1":
                pass
            elif key == "SCALE2":
                pass
            elif key == "SCALE3":
                pass
            elif key == "MTRIX1":
                pass
            elif key == "MTRIX2":
                pass
            elif key == "MTRIX3":
                pass
            elif key == "MODEL":
                pass
            elif key == "ATOM" or key == "HETATM":
                # serial = int(line[6:11])  # noqa: F841
                name = line[12:16].strip()
                # altloc = line[16]  # noqa: F841
                resname = line[17:20].strip()
                chainid = line[21]
                resseq = int(line[22:26])
                # icode = line[26]  # noqa: F841
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                tmp = line[54:60].strip()
                occupancy = 1.0 if tmp == "" else float(tmp)
                tmp = line[60:66].strip()
                tempfactor = 0.0 if tmp == "" else float(tmp)

                # Symbol maybe fully capitalized e.g. 'FE', so need to fix
                symbol = line[75:78].strip().capitalize()
                # tmp = line[78:80].strip()
                # charge = 0.0 if tmp == '' else float(tmp)  # noqa: F841

                names.append(name)
                symbols.append(symbol)
                xs.append(x)
                ys.append(y)
                zs.append(z)
                resnames.append(resname)
                chainids.append(chainid)
                resseqs.append(resseq)
                occupancies.append(occupancy)
                tempfactors.append(tempfactor)
            elif key == "ANISOU":
                pass
            elif key == "TER":
                pass
            elif key == "ENDMDL":
                pass
            elif key == "CONECT":
                if connections is None:
                    n_atoms = len(symbols)
                    connections = []
                    for i in range(n_atoms + 1):
                        connections.append(list())

                atom = int(line[6:11])
                # should use columns in case they run together.
                for i in range(11, 31, 5):
                    tmp = line[i : i + 5].strip()
                    if tmp == "":
                        break
                    connections[atom].append(int(tmp))
            elif key == "MASTER":
                pass
            elif key == "END":
                break
            else:
                raise RuntimeError("Illegal line in PDB file\n\t" + line)

        if "name" not in self.atoms:
            self.atoms.add_attribute("name", coltype="string")

        atom_id = self.atoms.append(symbol=symbols, name=names, x=xs, y=ys, z=zs)

        if "resname" in self.atoms:
            self.atoms["resname"][0:] = resnames
        else:
            counts = collections.Counter(resnames)
            if len(counts) > 1 or [*counts.keys()] != ["UNK"]:
                self.atoms.add_attribute("resname", coltype="string")
                self.atoms["resname"][0:] = resnames

        if "chainid" in self.atoms:
            self.atoms["chainid"][0:] = chainids
        else:
            counts = collections.Counter(chainids)
            if len(counts) > 1 or [*counts.keys()] != ["A"]:
                self.atoms.add_attribute("chainid", coltype="string")
                self.atoms["chainid"][0:] = chainids

        if "resseq" in self.atoms:
            self.atoms["resseq"][0:] = resseqs
        else:
            counts = collections.Counter(resseqs)
            if len(counts) > 1 or [*counts.keys()] != ["1"]:
                self.atoms.add_attribute("resseq", coltype="int", default=1)
                self.atoms["resseq"][0:] = resseqs

        if "occupancy" in self.atoms:
            self.atoms["occupancy"][0:] = occupancies
        else:
            counts = collections.Counter(occupancies)
            if len(counts) > 1 or [*counts.keys()] != [1.0]:
                self.atoms.add_attribute("occupancy", coltype="float", default=1.0)
                self.atoms["occupancy"][0:] = occupancies

        if "tempfactor" in self.atoms:
            self.atoms["tempfactor"][0:] = tempfactors
        else:
            counts = collections.Counter(tempfactors)
            if len(counts) > 1 or [*counts.keys()] != [0.0]:
                self.atoms.add_attribute("tempfactor", coltype="float", default=1.0)
                self.atoms["tempfactor"][0:] = tempfactors

        iatom = []
        jatom = []
        if connections is not None:
            for i in range(1, n_atoms + 1):
                for j in connections[i]:
                    if i not in connections[j]:
                        logger.warning(f"Bond {j}-{i} not found in PDB file")
                        # put in the bond since we won't see its partner!
                        iatom.append(atom_id[i - 1])
                        jatom.append(atom_id[j - 1])
                    elif i < j:
                        # put in only 1 of 2 equivalent bonds
                        iatom.append(atom_id[i - 1])
                        jatom.append(atom_id[j - 1])

            self.bonds.append(i=iatom, j=jatom)

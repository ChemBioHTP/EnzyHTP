"""Definition for the PDBLine class. This is a convenient representation of a PDB line that allows
for easy editing of raw files.

PDBLine objects should not be handled by non-developer users and are primarily intended as a lighter
weight representation of PDB objects when manipulating the files directly.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
# TODO(CJ): finish the documentation in here
from typing import List

import enzy_htp.core as core
from enzy_htp.core import file_system as fs
from enzy_htp.chemical import THREE_LETTER_AA_MAPPER


class PDBLine:
    """
    Represents a line in a PDB file. Stores the original line and provides a basic description
        of the data it holds. Has a number of getters for determining the type of line that
        it is.


    Attributes:
            line : The original line from the file, without newline.
                line_type : Type of PDB line/record.
                atom_id : One-indexed index of the atom.
                atom_name : Name of the atom as a string.
                resi_name : Name of the residue as a string.
                resi_id : One-indexed index of the atom.
                chain_id : Chain id, typically a single character/string of length 1.
                atom_x : x coordinate value with 3 points of decimal precision.
                atom_y : y coordinate value with 3 points of decimal precision.
                atom_z : z coordinate value with 3 points of decimal precision.
                                charge : The int() charge of the line.
    """

    def __init__(self, line: str):
        """
        Initilize with a specific line in the PDB file. Tries to parse out all information but stores the original str() as well.
        """
        self.resi_name = str
        self.line = line
        if not (self.is_ATOM() or self.is_HETATM()):
            return

        self.line_type = self.line[0:6].strip()  # atom
        self.atom_id = int(self.line[6:11])
        self.atom_name = self.line[12:16].strip()
        # residue
        self.resi_name = self.line[17:20].strip()
        self.resi_id = int(self.line[22:26])
        # chain
        self.chain_id = self.line[21:22]
        # coord
        self.atom_x = float(self.line[30:38])
        self.atom_y = float(self.line[38:46])
        self.atom_z = float(self.line[46:54])

    def is_REMARK(self) -> bool:
        """Checks if a line is a REMARK record."""
        return self.line.startswith("REMARK")

    def is_DBREF(self) -> bool:
        """Checks if a line is a DBREF record."""
        return self.line.startswith("DBREF")

    def is_JRNL(self) -> bool:
        """Checks if a line is a JRNL record."""
        return self.line.startswith("JRNL")

    def is_REVDAT(self) -> bool:
        """Checks if a line is a REVDAT record."""
        return self.line.startswith("REVDAT")

    def is_AUTHOR(self) -> bool:
        """Checks if a line is a AUTHOR record."""
        return self.line.startswith("AUTHOR")

    def is_EXPDTA(self) -> bool:
        """Checks if a line is a EXPDTA record."""
        return self.line.startswith("EXPDTA")

    def is_KEYWDS(self) -> bool:
        """Checks if a line is a KEYWDS record."""
        return self.line.startswith("KEYWDS")

    def is_SOURCE(self) -> bool:
        """Checks if a line is a SOURCE record."""
        return self.line.startswith("SOURCE")

    def build(self, ff: str = "AMBER") -> str:
        """Converts the PDBLine() object into a PDB line ready to be written to file. Is compatible with the specified force field (ff)."""
        # TODO(CJ): This function needs some work.
        # 1. Needs an F string instead of the weird concatentation
        # 2. throw some kind of error if the ff is not suported
        if not (self.is_ATOM() or self.is_HETATM()):
            return self.line
        self.get_alternate_location_indicator()
        self.get_charge()
        self.get_element()
        self.get_insert_code()
        self.get_occupancy()
        self.get_seg_id()
        self.get_temp_factor()

        if ff == "AMBER":
            # build an amber style line here
            l_type = "{:<6}".format(self.line_type)
            a_index = "{:>5d}".format(self.atom_id)

            if len(self.atom_name) > 3:
                a_name = "{:<4}".format(self.atom_name)
            else:
                a_name = "{:<3}".format(self.atom_name)
                a_name = " " + a_name
            r_name = f"{self.resi_name}"

            c_index = self.chain_id
            r_index = "{:>4d}".format(self.resi_id)
            x = "{:>8.3f}".format(self.atom_x)
            y = "{:>8.3f}".format(self.atom_y)
            z = "{:>8.3f}".format(self.atom_z)

            AL_id = "{:1}".format(self.AL_id)
            insert_code = "{:1}".format(self.insert_code)
            occupancy = "{:>6.2f}".format(self.occupancy)
            temp_factor = "{:>6.2f}".format(self.temp_factor)
            seg_id = "{:<4}".format(self.seg_id)
            element = "{:1}".format(self.element)
            charge = "{:2}".format(self.charge)

        # example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        return f"{l_type}{a_index} {a_name}{AL_id}{r_name} {c_index}{r_index}{insert_code}   {x}{y}{z}{occupancy}{temp_factor}       {seg_id}{element}{charge}"

    def is_COMPND(self) -> bool:
        """Checks if a line is a COMPND record."""
        return self.line.startswith("COMPND")

    def is_TITLE(self) -> bool:
        """Checks if a line is a TITLE record."""
        return self.line.startswith("TITLE")

    def is_HEADER(self) -> bool:
        """Checks if a line is a HEADER record."""
        return self.line.startswith("HEADER")

    def is_HETATM(self) -> bool:
        """Checks if a line is a HETATM record."""
        return self.line.startswith("HETATM")

    def is_ATOM(self) -> bool:
        """Checks if a line is a ATOM record."""
        return self.line.startswith("ATOM")

    def is_TER(self) -> bool:
        """Checks if a line is a "TER" terminating code."""
        return self.line.startswith("TER")

    def is_END(self) -> bool:
        """Checks if a line is an "END" end of file code."""
        return self.line[:3] == "END"

    def is_CRYST1(self) -> bool:
        """Checks if a line starts with a "CRYST1" code."""
        return self.line[:6] == "CRYST1"

    def is_residue_line(self) -> bool:
        """Checks if a line is part of a residue."""
        return self.resi_name in THREE_LETTER_AA_MAPPER

    # misc
    def get_alternate_location_indicator(self):
        self.AL_id = self.line[16:17].strip()
        return self.AL_id

    def get_insert_code(self):
        self.insert_code = self.line[26:27].strip()
        return self.insert_code

    def get_occupancy(self):
        self.occupancy = self.line[54:60].strip()
        if self.occupancy == "":
            self.occupancy = 1.0
        else:
            self.occupancy = float(self.occupancy)
        return self.occupancy

    def get_temp_factor(self):
        self.temp_factor = self.line[60:66].strip()
        if self.temp_factor == "":
            self.temp_factor = 0.0
        else:
            self.temp_factor = float(self.temp_factor)
        return self.temp_factor

    def get_seg_id(self):
        self.seg_id = self.line[72:76].strip()
        return self.seg_id

    def get_element(self) -> str:
        self.element = self.line[76:78].strip()
        # try to get if not exist
        # if self.element == '':
        #     if self.atom_name in Resi_Ele_map['Amber']:
        #         self.element = Resi_Ele_map['Amber'][self.atom_name]
        return self.element

    def is_water(self) -> bool:
        """Checks if a line is an alias to water."""
        return self.resi_name in {"Na+", "Cl-", "WAT", "HOH"}

    def get_charge(self) -> str:
        """Finds and returns the charge of line as a string."""
        self.charge = self.line[78:80].strip()
        return self.charge

    def __str__(self) -> str:
        return self.line

    def __repr__(self) -> str:
        return self.line


def read_pdb_lines(fname: str) -> List[PDBLine]:
    """Generates a list() of PDBLine objexts from a given PDB file"""
    # TODO(CJ): there is a chance that we only want HETATM/ATOM/TER/END... will figure out if this is the case later
    ending = fs.get_file_ext(fname).lower()
    if not ending == ".pdb":
        raise core.UnsupportedFileType(
            f"read_pdb_lines() requires a .pdb file. Could not read '{fname}'"
        )
    non_empty = list(filter(len, fs.lines_from_file(fname)))
    return list(map(PDBLine, non_empty))

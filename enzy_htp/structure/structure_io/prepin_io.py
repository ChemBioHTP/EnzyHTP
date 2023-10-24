"""Generation/construction of Structure() objects from .prepin files and exporting these objects to this file format. 
Definition of .prepin format (https://ambermd.org/doc/prep.html). All parsing is done within enzy_htp using 
this parser only. The PrepinParser has no private data and serves as a namespace for .prepin I/O conversion functions.

Author: Qianzhen Shao <shaoqz@icloud.com>
Date: 2023-10-17
"""
from typing import Dict, List

from ._interface import StructureParserInterface
from ..structure import Structure, convert_res_to_structure
from ..atom import Atom
from ..residue import Residue

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import FileFormatError

class PrepinParser(StructureParserInterface):
    """the parser for AmberMD prepin files"""

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """pass"""
        pass

    @classmethod
    def get_structure(cls, path: str) -> Structure:
        """Converting a .prepin file (as its path) into the Structure()
        Arg:
            path:
                the file path of the PDB file
        Return:
            Structure()
        """
        prepin_data = cls._parse_prepin_file(path)

        # support check
        if prepin_data["IFIXC"] == "CHANGE":
            _LOGGER.error("Does not support parse CHANGE prepin file to Structure() yet. "
                          f"Please post an issue if really needed. ({path})")
            raise FileFormatError

        # build Structure()
        atoms = cls._build_atoms(prepin_data)
        cls._connect_atoms(atoms, prepin_data)
        res = cls._build_residue(atoms, prepin_data)
        return convert_res_to_structure(res)

    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        """convert a Structure() to .prepin file content. Only 1 residue unit is allowed in the stru"""
        pass # TODO

    @classmethod
    def _parse_prepin_file(cls, path: str) -> Dict:
        """parse prepin file to a data dictionary. The name of keys are based on
        'Input  description' section of https://ambermd.org/doc/prep.html with
        slight modification described in docstrings of called function below"""
        if not fs.has_content( path ):
            _LOGGER.error(f"The supplied file {path} does not exist or is empty. Exiting...")
            raise ValueError
                
        result = {}
        lines: List[str] = fs.lines_from_file( path )

        # Section -1-
        result.update(cls._parse_database_control(lines[0]))
        # Section -2-
        result.update(cls._parse_namdbf(lines[1]))
        # Section -3-
        result.update(cls._parse_title(lines[2]))
        # Section -4-
        result.update(cls._parse_namf(lines[3]))
        # Section -5-
        result.update(cls._parse_resname(lines[4]))
        # Section -6-
        result.update(cls._parse_icontrol(lines[5]))
        # Section -7-
        result.update(cls._parse_cut(lines[6]))
        # Section -8-
        coord_end = cls._deduce_coord_end(lines)
        result.update(cls._parse_coord(lines[7:coord_end]))
        # Section -9-
        result.update(cls._parse_iopr(lines[coord_end:]))

        return result

    @classmethod
    def _parse_database_control(cls, line: str) -> Dict:
        """Fortran format: FORMAT(3I). Each means:
        IDBGEN:
            = 0  No database generation.  Output will be individual files.
                This is the standard procedure if you want to create a
                single small molecule.
            = 1  A new data base will be generated or the existing database
                will be appended.
        IREST       Flag for the type of generation (assuming IDBGEN = 1)
            = 0  New data base
            = 1  Appending an existing data base

        ITYPF       Force field type code (used in LINK stage)
            Ignored if IDBGEN = 0   The following codes are used in
            the standard database:
            = 1  United atom model
            = 2  All atom model
            = 100  United atom charged N-terminal amino acid residues
            = 101  United atom charged C-terminal amino acid residues
            = 200  All atom charged N-terminal amino acid residues
            = 201  All atom charged C-terminal amino acid residues"""
        line_parts = line.strip().split()
        if len(line_parts) != 3:
            _LOGGER.error(f"Invalid prepin format! '{line}'")
            raise FileFormatError
        idbgen, irest, itypf = line_parts

        return {
            "IDBGEN" : idbgen,
            "IREST" : irest,
            "ITYPF" : itypf,}

    @classmethod
    def _parse_namdbf(cls, line: str) -> Dict:
        """Fortran format: FORMAT(A80).
        NAMDBF      Name of the data base file (maximum 80 characters)
            if NOT data base generation leave a BLANK CARD"""
        return {"NAMDBF" : line.strip()}

    @classmethod
    def _parse_title(cls, line: str) -> Dict:
        """Fortran format: FORMAT(20A4) #seems not?
        TITLE      Descriptive header for the residue"""
        return {"TITLE" : line.strip()}  

    @classmethod
    def _parse_namf(cls, line: str) -> Dict:
        """Fortran format: FORMAT(A80)
        NAMF       Name of the output file if an individual residue file is
            being generated.  If database is being generated or
            appended this card IS read but ignored."""
        return {"NAMF" : line.strip()}

    @classmethod
    def _parse_resname(cls, line: str) -> Dict:
        """Fortran format: FORMAT(2A,I)
        NAMRES     A unique name for the residue of maximum 4 characters
        INTX       Flag for the type of coordinates to be saved for the
            LINK module
            'INT'  internal coordinates will be output (preferable)
            'XYZ'  cartesian coordinates will be output

        KFORM      Format of output for individual residue files
            = 0  formatted output (recommended for debugging)
            = 1  binary output"""
        line_parts = line.strip().split()
        if len(line_parts) != 3:
            _LOGGER.error(f"Invalid prepin format! '{line}'")
            raise FileFormatError

        namres, intx, kform = line_parts

        return {
            "NAMRES" : namres,
            "INTX" : intx,
            "KFORM" : kform,}

    @classmethod
    def _parse_icontrol(cls, line: str) -> Dict:
        """Fortran format: FORMAT(4A)
        IFIXC      Flag for the type of input geometry of the residue(s)

            'CORRECT' The geometry is input as internal coordinates with
                        correct order according to the tree structure.
                        NOTE: the tree structure types ('M', 'S', etc) and order
                        must be defined correctly: NA(I), NB(I), and NC(I) on card
                        8 are always ignored.
            'CHANGE'  It is input as cartesian coordinates or part cartesian
                        and part internal.  Cartesians should precede internals
                        to ensure that the resulting coordinates are correct.
                        Coordinates need not be in correct order, since each
                        is labeled with its atom number. NOTE: NA(I), NB(I), and
                        NC(I) on card 8 must be omitted for cartesian coordinates
                        with this option.

        IOMIT      Flag for the omission of dummy atoms

            'OMIT'    dummy atoms will be deleted after generating all the
                        information (this is used for all but the first residue
                        in the system)
            'NOMIT'   they will not be deleted (dummy atoms are retained for
                        the first residue of the system.  others are omitted)

        ISYMDU     Symbol for the dummy atoms.  The symbol must be
             be unique.  It is preferable to use 'DU' for it

        IPOS       Flag for the position of dummy atoms to be deleted

            'ALL'     all the dummy atoms will be deleted
            'BEG'     only the beginning dummy atoms will be deleted"""
        line_parts = line.strip().split()
        if len(line_parts) != 4:
            _LOGGER.error(f"Invalid prepin format! '{line}'")
            raise FileFormatError

        ifixc, iomit, isymdu, ipos = line_parts

        return {
            "IFIXC" : ifixc,
            "IOMIT" : iomit,
            "ISYMDU" : isymdu,
            "IPOS" : ipos, }

    @classmethod
    def _parse_cut(cls, line: str) -> Dict:
        """Fortran format: FORMAT(F)
        CUT  The cutoff distance for loop closing bonds which
             cannot be defined by the tree structure.  Any pair of
             atoms within this distance is assumed to be bonded.
             We recommend that CUT be set to 0.0 and explicit loop
             closing bonds be defined below"""
        return {"CUT" : line.strip()}

    @classmethod
    def _deduce_coord_end(cls, lines: List[str]) -> int:
        """determine the ending line index of the coordinate section.
        This section is terminated by one or two blank lines.
        Return:
            int, the index of the blank line at the end of the coord
            section."""
        for i, line in enumerate(lines):
            if i > 6:
                if not line.strip():
                    return i

    @classmethod
    def _parse_coord(cls, lines: List[str]) -> Dict:
        """coordinate_section: List[Dict], each list element comes from each line.
        For each line:
        Fortran format: FORMAT(I,3A,3I,4F)
            id
                The actual number of the atom in the tree.

                If IFIXC .eq. 'CHANGE' then this number is important
                since the corresponding coordinates are stored at that
                location.  If IFIXC .eq. 'CORRECT' then atoms are in
                the correct order according to the tree structure.

            NOTE:  PREP always expects three dummy atoms for the beginning.

            name
                A unique atom name for the atom I. If coordinates are
                read in at the EDIT stage, this name will be used for
                matching atoms.  Maximum 4 characters.

            atomtype
                A symbol for the atom I which defines its force field
                atom type and is used in the module PARM for assigning
                the force field parameters.

            ITREE
                The topological type (tree symbol) for atom I
                (M, S, B, E, or 3)

            NA(I)      The atom number to which atom I is connected.
                Read but ignored for internal coordinates; If cartesian
                coordinates are used, this must be omitted.

            NB(I)      The atom number to which atom I makes an angle along
                with NA(I).
                Read but ignored for internal coordinates; If cartesian
                coordinates are used, this must be omitted.

            NC(I)      The atom number to which atom I makes a dihedral along
                with NA(I) and NB(I).
                Read but ignored for internal coordinates; If cartesian
                coordinates are used, this must be omitted.

            R(I)       If IFIXC .eq. 'CORRECT' then this is the bond length
                between atoms I and NA(I)
                If IFIXC .eq. 'CHANGE' then this is the X coordinate
                of atom I

            THETA(I)   If IFIXC .eq. 'CORRECT' then it is the bond angle
                between atom NB(I), NA(I) and I
                If IFIXC .eq. 'CHANGE' then it is the Y coordinate of
                atom I

            PHI(I)     If IFIXC .eq. 'CORRECT' then it is the dihedral angle
                between NC(I), NB(I), NA(I) and I
                If IFIXC .eq. 'CHANGE' then it is the Z coordinate of
                atom I

            CHRG(I)    The partial atomic charge on atom I

            This section is terminated by one BLANK CARD if IFIXC = 'CORRECT'.
            This section is terminated by TWO BLANK CARDS if IFIXC = 'CHANGE'
            # seems all use  2 empty lines now"""
        result = []
        for line in lines:
            line_parts = line.strip().split()
            if len(line_parts) != 11:
                _LOGGER.error(f"Invalid prepin format! '{line}'")
                raise FileFormatError

            result.append({
                "id" : line_parts[0],
                "name" : line_parts[1],
                "atomtype" : line_parts[2],
                "ITREE" : line_parts[3],
                "NA" : line_parts[4],
                "NB" : line_parts[5],
                "NC" : line_parts[6],
                "R" : line_parts[7],
                "THETA" : line_parts[8],
                "PHI" : line_parts[9],
                "CHRG" : line_parts[10],
            })

        return {"coordinate_section" : result}

    @classmethod
    def _parse_iopr(cls, lines: List[str]) -> Dict:
        """Fortran format: FORMAT(A4)
        IOPR
            Flag to read additional information about the residue.
            There are four options available.  The order in which
            they are specified is not important.  Format is keyword
            on its own line, followed by data on succeeding lines,
            terminated by a BLANK CARD."""
        result = {}
        # divide sections
        loop_lines = []
        charge_lines = []
        improper_lines = []
        for i, line in enumerate(lines):
            line = line.strip()
            if line == "LOOP":
                for belowline in lines[i+1:]:
                    if belowline.strip() == "":
                        break
                    loop_lines.append(belowline)
            if line == "IMPROPER":
                for belowline in lines[i+1:]:
                    if belowline.strip() == "":
                        break
                    improper_lines.append(belowline)
            if line == "CHARGE":
                for belowline in lines[i+1:]:
                    if belowline.strip() == "":
                        break
                    charge_lines.append(belowline)
            if line == "DONE":
                break
        if loop_lines:
            result["LOOP"] = cls._parse_loop(loop_lines)
        if improper_lines:
            result["IMPROPER"] = cls._parse_improper(improper_lines)
        if charge_lines:
            result["CHARGE"] = cls._parse_charge(charge_lines)
        return result

    @classmethod
    def _parse_loop(cls, loop_lines: List[str]) -> List:
        """'LOOP'
        Fortran format: FORMAT(2A)
            Control to read explicit loop closing bonds (in
            addition to the loops generated based on the cutoff
            criterion).  If this option is used it is preferable
            to set the cutoff criterion to zero.  The loop closing
            atoms are read in format(2A) as their atom (IGRAPH) names.
            A BLANK CARD terminates this section.
        Return:
            [[a1_name, a2_name], ...]"""
        result = []
        for line in loop_lines:
            line_parts = line.strip().split()
            if len(line_parts) != 2:
                _LOGGER.error(f"Invalid prepin format! '{line}'")
                raise FileFormatError

            result.append(line_parts)

        return result

    @classmethod
    def _parse_improper(cls, improper_lines: List[str]) -> Dict:
        """'IMPROPER'
        Fortran format: FORMAT(4A)
            Control for reading the improper torsion angles.  A
            proper torsion I - J - K - L has I bonded to J bonded
            to K bonded to L.  An IMPROPER torsion is any torsion in
            which this is not the case.  Improper torsions are used to
            keep the asymmetric centers from racemizing in the united
            atom model where all the C-H hydrogens are omitted.  They
            can also be used to enforce planarity.  The normal case is:

                               J
                               |
                               K
                              / \
                             I   L

                        Improper I-J-K-L

            where the central atom (K) is the third atom in the improper
            and the order of the other three is determined alphabetically
            by atom type and if types are the same by atom number.
            The improper torsions should be defined in such a way that
            the proper torsions are not duplicated.  The atoms making the
            improper torsions are read as their atom (IGRAPH) names.
            '-M' can be used in place of an atom name to indicate the
            last main chain atom in the previous residue, and '+M' for
            the first main chain atom in the next residue. NOTE: -M and
            +M cannot be used in the 4th position ('L') owing to internal
            data representation limitations.  A BLANK CARD terminates
            this section.
        # seems not useful in EnzyHTP"""
        result = []
        for line in improper_lines:
            line_parts = line.strip().split()
            if len(line_parts) != 4:
                _LOGGER.error(f"Invalid prepin format! '{line}'")
                raise FileFormatError

            result.append(line_parts)

        return result

    @classmethod
    def _parse_charge(cls, charge_lines: List[str]) -> List:
        """'CHARGE'
        Fortran format: FORMAT(5F)
            Control to read additional partial atomic charges.
            These will override charges specified above in section 8.
            The charges are read in format(5F) for the non-dummy
            atoms.  A BLANK CARD terminates this section.   It is
            less error-prone to specify charges as in section 8.
        # Are they for all atoms?"""
        _LOGGER.error("Add support for this part in PrepinParser._parse_charge!")
        raise Exception

    @classmethod
    def _build_atom(cls, prepin_data: Dict) -> List[Atom]:
        """build a list of Atom()s from {prepin_data}"""
        labels:List[str]="x_coord y_coord z_coord atom_name charge atom_type".split()
        temp = dict()
        for ll in labels:
            temp[ll] = aa[ll]
        # TODO

    @classmethod 
    def _connect_atoms(cls, atoms: List[Atom], prepin_data: Dict):
        """connect {atoms} using information from {prepin_data}.
        connectivities are written to {atom.connect} of each {atom}"""

    @classmethod 
    def _build_residue(cls, atoms: List[Atom], prepin_data: Dict) -> Residue:
        """build Residue() based on {atoms} and {prepin_data}"""

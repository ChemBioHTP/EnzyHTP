"""Stores mappers and testing functions for residue information. Much of this functionality is later used in 
NOTE: using PDB style atom names as official atom naming scheme in EnzyHTP. Make atom name convertor for other
    naming scheme.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from typing import List, Set, Dict, Union
from ..core import InvalidResidueCode, _LOGGER

# yapf: disable

AA_LIST: List[str] = ["R","H","K","D","E","S","T","N","Q","C","G","P","A","V","I","L","M","F","Y","W","U"]
"""Capitalized list of all one-letter amino acid names."""

THREE_LETTER_AA_MAPPER: Dict[str, str] = {
        "ALA": "A",
        "ARG": "R",
        "ASH": "N",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "CYM": "C",
        "CYX": "C",
        "GLN": "Q",
        "GLH": "E",
        "GLU": "E",
        "GLY": "G",
        "HID": "H",
        "HIE": "H",
        "HIP": "H",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYN": "K",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SEC": "U",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V"
} #TODO(shaoqz): @imp2 add related to canonical in the name
"""Contains mapping of all amino acids codes, with key value pairs of (three letter code, one letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_one_letter()"""

ONE_LETTER_AA_MAPPER: Dict[str, str] = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "U": "SEC",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR"
}
"""Contains mapping of all amino acids codes, with key value pairs of (one letter code, three letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_three_letter()"""

ONE_LETTER_CAA_MAPPER: Dict[str, str] = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR"
}
"""Contains mapping of all canonical amino acids codes, with key value pairs of (one letter code, three letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_one_letter()"""

THREE_TO_THREE_LETTER_CAA_MAPPER: Dict[str, str] = {
        "ALA": "ALA",
        "ARG": "ARG",
        "ASH": "ASN",
        "ASN": "ASN",
        "ASP": "ASP",
        "CYS": "CYS",
        "CYM": "CYS",
        "CYX": "CYS",
        "GLN": "GLN",
        "GLH": "GLU",
        "GLU": "GLU",
        "GLY": "GLY",
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "HIS": "HIS",
        "ILE": "ILE",
        "LEU": "LEU",
        "LYN": "LYS",
        "LYS": "LYS",
        "MET": "MET",
        "PHE": "PHE",
        "PRO": "PRO",
        "SEC": "SEC",
        "SER": "SER",
        "THR": "THR",
        "TRP": "TRP",
        "TYR": "TYR",
        "VAL": "VAL"
}
"""Contains mapping of all canonical amino acids codes, with key value pairs of
(three letter code, canonicalized three letter code). Should NOT be called directly
for code conversion. Instead used enzy_htp.chemical.residue.convert_to_canoncial_three_letter()"""

RESIDUE_ELEMENT_MAP: Dict[str, Dict[str, str]] = {
                "Amber": {    "C" : "C",
                              "CA" : "C",
                              "CB" : "C",
                              "CD" : "C",
                              "CD1" : "C",
                              "CD2" : "C",
                              "CE" : "C",
                              "CE1" : "C",
                              "CE2" : "C",
                              "CE3" : "C",
                              "CG" : "C",
                              "CG1" : "C",
                              "CG2" : "C",
                              "CH2" : "C",
                              "CZ" : "C",
                              "CZ2" : "C",
                              "CZ3" : "C",
                              "H" : "H",
                              "H1" : "H",
                              "H2" : "H",
                              "H3" : "H",
                              "HA" : "H",
                              "HA2" : "H",
                              "HA3" : "H",
                              "HB" : "H",
                              "HB1" : "H",
                              "HB2" : "H",
                              "HB3" : "H",
                              "HD1" : "H",
                              "HD11" : "H",
                              "HD12" : "H",
                              "HD13" : "H",
                              "HD2" : "H",
                              "HD21" : "H",
                              "HD22" : "H",
                              "HD23" : "H",
                              "HD3" : "H",
                              "HE" : "H",
                              "HE1" : "H",
                              "HE2" : "H",
                              "HE21" : "H",
                              "HE22" : "H",
                              "HE3" : "H",
                              "HG" : "H",
                              "HG1" : "H",
                              "HG11" : "H",
                              "HG12" : "H",
                              "HG13" : "H",
                              "HG2" : "H",
                              "HG21" : "H",
                              "HG22" : "H",
                              "HG23" : "H",
                              "HG3" : "H",
                              "HH" : "H",
                              "HH11" : "H",
                              "HH12" : "H",
                              "HH2" : "H",
                              "HH21" : "H",
                              "HH22" : "H",
                              "HZ" : "H",
                              "HZ1" : "H",
                              "HZ2" : "H",
                              "HZ3" : "H",
                              "N" : "N",
                              "ND1" : "N",
                              "ND2" : "N",
                              "NE" : "N",
                              "NE1" : "N",
                              "NE2" : "N",
                              "NH1" : "N",
                              "NH2" : "N",
                              "NZ" : "N",
                              "O" : "O",
                              "OD1" : "O",
                              "OD2" : "O",
                              "OE1" : "O",
                              "OE2" : "O",
                              "OG" : "O",
                              "OG1" : "O",
                              "OH" : "O",
                              "OXT" : "O",
                              "SD" : "S",
                              "SG" : "S",
                              "LI" : "Li",
                              "NA" : "Na",
                              "K"  : "K" ,
                              "RB" : "Rb",
                              "CS" : "Cs",
                              "MG" : "Mg",
                              "TL" : "Tl",
                              "CU" : "Cu",
                              "AG" : "Ag",
                              "BE" : "Be",
                              "NI" : "Ni",
                              "PT" : "Pt",
                              "ZN" : "Zn",
                              "CO" : "Co",
                              "PD" : "Pd",
                              "CR" : "Cr",
                              "FE" : "Fe",
                              "V"  : "V" ,
                              "MN" : "Mn",
                              "YB" : "Yb",
                              "SN" : "Sn",
                              "PB" : "Pb",
                              "EU" : "Eu",
                              "SR" : "Sr",
                              "SM" : "Sm",
                              "BA" : "Ba",
                              "RA" : "Ra",
                              "AL" : "Al",
                              "IN" : "In",
                              "Y"  : "Y" ,
                              "LA" : "La",
                              "PR" : "Pr",
                              "ND" : "Nd",
                              "GD" : "Gd",
                              "TB" : "Tb",
                              "DY" : "Dy",
                              "ER" : "Er",
                              "TM" : "Tm",
                              "LU" : "Lu",
                              "HF" : "Hf",
                              "ZR" : "Zr",
                              "U"  : "U" ,
                              "PU" : "Pu",
                              "TH" : "Th" }
}
"""Mapper that shows element names for alternative names of atoms. Key value structure is (force field name, (atom mapper)), where atom mapper is an additional
mapper that maps altered atom names to base atom names. Currently only defined for "Amber".

Example usage:

>>> initial_atom_name = "CA"
>>> RESIDUE_ELEMENT_MAP["Amber"][initial_atom_name] TODO also remove those ff keys
"C"

"""

RESIDUE_CONNECTIVITY_MAP: Dict[str, Dict[str, List[str]]] = {
                  "ALA": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB1", "HB2", "HB3"], "HB1": ["CB"], "HB2": ["CB"], "HB3": ["CB"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "ARG": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "NE"], "HD2": ["CD"], "HD3": ["CD"], "NE": ["CD", "HE", "CZ"], "HE": ["NE"], "CZ": ["NE", "NH1", "NH2"], "NH1": ["CZ", "HH11", "HH12"], "HH11": ["NH1"], "HH12": ["NH1"], "NH2": ["CZ", "HH21", "HH22"], "HH21": ["NH2"], "HH22": ["NH2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "ASH": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "OD2"], "OD1": ["CG"], "OD2": ["CG", "HD2"], "HD2": ["OD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "ASN": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "ND2"], "OD1": ["CG"], "ND2": ["CG", "HD21", "HD22"], "HD21": ["ND2"], "HD22": ["ND2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "ASP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "OD2"], "OD1": ["CG"], "OD2": ["CG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "CYM": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB3", "HB2", "SG"], "HB3": ["CB"], "HB2": ["CB"], "SG": ["CB"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "CYS": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB", "HG"], "HG": ["SG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "CYX": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "GLH": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "OE2"], "OE1": ["CD"], "OE2": ["CD", "HE2"], "HE2": ["OE2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "GLN": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "NE2"], "OE1": ["CD"], "NE2": ["CD", "HE21", "HE22"], "HE21": ["NE2"], "HE22": ["NE2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "GLU": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "OE2"], "OE1": ["CD"], "OE2": ["CD"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "GLY": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA2", "HA3", "C"], "HA2": ["CA"], "HA3": ["CA"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "HID": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "CD2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "HIE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "CE1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "HIP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "HYP": {"N": ["CD", "CA", "-1C"], "CD": ["N", "HD22", "HD23", "CG"], "HD22": ["CD"], "HD23": ["CD"], "CG": ["CD", "HG", "OD1", "CB"], "HG": ["CG"], "OD1": ["CG", "HD1"], "HD1": ["OD1"], "CB": ["CG", "HB2", "HB3", "CA"], "HB2": ["CB"], "HB3": ["CB"], "CA": ["N", "CB", "HA", "C"], "HA": ["CA"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "ILE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "CG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "CG1": ["CB", "HG12", "HG13", "CD1"], "HG12": ["CG1"], "HG13": ["CG1"], "CD1": ["CG1", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "LEU": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG", "CD1", "CD2"], "HG": ["CG"], "CD1": ["CG", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "CD2": ["CG", "HD21", "HD22", "HD23"], "HD21": ["CD2"], "HD22": ["CD2"], "HD23": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "LYN": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "CE"], "HD2": ["CD"], "HD3": ["CD"], "CE": ["CD", "HE2", "HE3", "NZ"], "HE2": ["CE"], "HE3": ["CE"], "NZ": ["CE", "HZ2", "HZ3"], "HZ2": ["NZ"], "HZ3": ["NZ"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "LYS": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "CE"], "HD2": ["CD"], "HD3": ["CD"], "CE": ["CD", "HE2", "HE3", "NZ"], "HE2": ["CE"], "HE3": ["CE"], "NZ": ["CE", "HZ1", "HZ2", "HZ3"], "HZ1": ["NZ"], "HZ2": ["NZ"], "HZ3": ["NZ"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "MET": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "SD"], "HG2": ["CG"], "HG3": ["CG"], "SD": ["CG", "CE"], "CE": ["SD", "HE1", "HE2", "HE3"], "HE1": ["CE"], "HE2": ["CE"], "HE3": ["CE"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "PHE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "HZ", "CE2"], "HZ": ["CZ"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "PRO": {"N": ["CD", "CA", "-1C"], "CD": ["N", "HD2", "HD3", "CG"], "HD2": ["CD"], "HD3": ["CD"], "CG": ["CD", "HG2", "HG3", "CB"], "HG2": ["CG"], "HG3": ["CG"], "CB": ["CG", "HB2", "HB3", "CA"], "HB2": ["CB"], "HB3": ["CB"], "CA": ["N", "CB", "HA", "C"], "HA": ["CA"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "SER": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "OG"], "HB2": ["CB"], "HB3": ["CB"], "OG": ["CB", "HG"], "HG": ["OG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "THR": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "OG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "OG1": ["CB", "HG1"], "HG1": ["OG1"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "TRP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "NE1"], "HD1": ["CD1"], "NE1": ["CD1", "HE1", "CE2"], "HE1": ["NE1"], "CE2": ["NE1", "CZ2", "CD2"], "CZ2": ["CE2", "HZ2", "CH2"], "HZ2": ["CZ2"], "CH2": ["CZ2", "HH2", "CZ3"], "HH2": ["CH2"], "CZ3": ["CH2", "HZ3", "CE3"], "HZ3": ["CZ3"], "CE3": ["CZ3", "HE3", "CD2"], "HE3": ["CE3"], "CD2": ["CG", "CE2", "CE3"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "TYR": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "OH", "CE2"], "OH": ["CZ", "HH"], "HH": ["OH"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "VAL": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG1", "CG2"], "HB": ["CB"], "CG1": ["CB", "HG11", "HG12", "HG13"], "HG11": ["CG1"], "HG12": ["CG1"], "HG13": ["CG1"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                  "WAT": {"O": ["H1", "H2"], "H1":["O"], "H2": ["O"]}
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms)."""

RESIDUE_CONNECTIVITY_MAP_CTERMINAL: Dict[str, Dict[str, List[str]]] = {
                        "ALA": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB1", "HB2", "HB3"], "HB1": ["CB"], "HB2": ["CB"], "HB3": ["CB"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "ARG": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "NE"], "HD2": ["CD"], "HD3": ["CD"], "NE": ["CD", "HE", "CZ"], "HE": ["NE"], "CZ": ["NE", "NH1", "NH2"], "NH1": ["CZ", "HH11", "HH12"], "HH11": ["NH1"], "HH12": ["NH1"], "NH2": ["CZ", "HH21", "HH22"], "HH21": ["NH2"], "HH22": ["NH2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "ASN": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "ND2"], "OD1": ["CG"], "ND2": ["CG", "HD21", "HD22"], "HD21": ["ND2"], "HD22": ["ND2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "ASP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "OD2"], "OD1": ["CG"], "OD2": ["CG"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "CYS": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB", "HG"], "HG": ["SG"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "CYX": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "GLN": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "NE2"], "OE1": ["CD"], "NE2": ["CD", "HE21", "HE22"], "HE21": ["NE2"], "HE22": ["NE2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "GLU": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "OE2"], "OE1": ["CD"], "OE2": ["CD"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "GLY": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA2", "HA3", "C"], "HA2": ["CA"], "HA3": ["CA"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "HID": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "CD2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "HIE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "CE1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "HIP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "HYP": {"N": ["CD", "CA", "-1C"], "CD": ["N", "HD22", "HD23", "CG"], "HD22": ["CD"], "HD23": ["CD"], "CG": ["CD", "HG", "OD1", "CB"], "HG": ["CG"], "OD1": ["CG", "HD1"], "HD1": ["OD1"], "CB": ["CG", "HB2", "HB3", "CA"], "HB2": ["CB"], "HB3": ["CB"], "CA": ["N", "CB", "HA", "C"], "HA": ["CA"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "ILE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "CG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "CG1": ["CB", "HG12", "HG13", "CD1"], "HG12": ["CG1"], "HG13": ["CG1"], "CD1": ["CG1", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "LEU": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG", "CD1", "CD2"], "HG": ["CG"], "CD1": ["CG", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "CD2": ["CG", "HD21", "HD22", "HD23"], "HD21": ["CD2"], "HD22": ["CD2"], "HD23": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "LYS": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "CE"], "HD2": ["CD"], "HD3": ["CD"], "CE": ["CD", "HE2", "HE3", "NZ"], "HE2": ["CE"], "HE3": ["CE"], "NZ": ["CE", "HZ1", "HZ2", "HZ3"], "HZ1": ["NZ"], "HZ2": ["NZ"], "HZ3": ["NZ"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "MET": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "SD"], "HG2": ["CG"], "HG3": ["CG"], "SD": ["CG", "CE"], "CE": ["SD", "HE1", "HE2", "HE3"], "HE1": ["CE"], "HE2": ["CE"], "HE3": ["CE"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "PHE": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "HZ", "CE2"], "HZ": ["CZ"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "PRO": {"N": ["CD", "CA", "-1C"], "CD": ["N", "HD2", "HD3", "CG"], "HD2": ["CD"], "HD3": ["CD"], "CG": ["CD", "HG2", "HG3", "CB"], "HG2": ["CG"], "HG3": ["CG"], "CB": ["CG", "HB2", "HB3", "CA"], "HB2": ["CB"], "HB3": ["CB"], "CA": ["N", "CB", "HA", "C"], "HA": ["CA"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "SER": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "OG"], "HB2": ["CB"], "HB3": ["CB"], "OG": ["CB", "HG"], "HG": ["OG"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "THR": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "OG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "OG1": ["CB", "HG1"], "HG1": ["OG1"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "TRP": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "NE1"], "HD1": ["CD1"], "NE1": ["CD1", "HE1", "CE2"], "HE1": ["NE1"], "CE2": ["NE1", "CZ2", "CD2"], "CZ2": ["CE2", "HZ2", "CH2"], "HZ2": ["CZ2"], "CH2": ["CZ2", "HH2", "CZ3"], "HH2": ["CH2"], "CZ3": ["CH2", "HZ3", "CE3"], "HZ3": ["CZ3"], "CE3": ["CZ3", "HE3", "CD2"], "HE3": ["CE3"], "CD2": ["CG", "CE2", "CE3"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "TYR": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "OH", "CE2"], "OH": ["CZ", "HH"], "HH": ["OH"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "VAL": {"N": ["H", "CA", "-1C"], "H": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG1", "CG2"], "HB": ["CB"], "CG1": ["CB", "HG11", "HG12", "HG13"], "HG11": ["CG1"], "HG12": ["CG1"], "HG13": ["CG1"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "C": ["CA", "O", "OXT"], "O": ["C"], "OXT": ["C"]},
                        "NHE": {"N": ["HN1", "HN2", "-1C"], "HN1": ["N"], "HN2": ["N"]},
                        "NME": {"N": ["H", "CH3", "-1C"], "H": ["N"], "CH3": ["N", "HH31", "HH32", "HH33"], "HH31": ["CH3"], "HH32": ["CH3"], "HH33": ["CH3"]}
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the C-terminal version of each residue."""

RESIDUE_CTERMINAL_ATOM_LIST: Dict[str, List[str]] = {
                        "ALA": ["N", "H", "CA", "HA", "CB", "HB1", "HB2", "HB3", "C", "O", "OXT"],
                        "ARG": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22", "C", "O", "OXT"],
                        "ASN": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22", "C", "O", "OXT"],
                        "ASP": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "OD1", "OD2", "C", "O", "OXT"],
                        "CYS": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "SG", "HG", "C", "O", "OXT"],
                        "CYX": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "SG", "C", "O", "OXT"],
                        "GLN": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22", "C", "O", "OXT"],
                        "GLU": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2", "C", "O", "OXT"],
                        "GLY": ["N", "H", "CA", "HA2", "HA3", "C", "O", "OXT"],
                        "HID": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "CD2", "HD2", "C", "O", "OXT"],
                        "HIE": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "ND1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "O", "OXT"],
                        "HIP": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "O", "OXT"],
                        "HYP": ["N", "CD", "HD22", "HD23", "CG", "HG", "OD1", "HD1", "CB", "HB2", "HB3", "CA", "HA", "C", "O", "OXT"],
                        "ILE": ["N", "H", "CA", "HA", "CB", "HB", "CG2", "HG21", "HG22", "HG23", "CG1", "HG12", "HG13", "CD1", "HD11", "HD12", "HD13", "C", "O", "OXT"],
                        "LEU": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23", "C", "O", "OXT"],
                        "LYS": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3", "NZ", "HZ1", "HZ2", "HZ3", "C", "O", "OXT"],
                        "MET": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "SD", "CE", "HE1", "HE2", "HE3", "C", "O", "OXT"],
                        "PHE": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "HZ", "CE2", "HE2", "CD2", "HD2", "C", "O", "OXT"],
                        "PRO": ["N", "CD", "HD2", "HD3", "CG", "HG2", "HG3", "CB", "HB2", "HB3", "CA", "HA", "C", "O", "OXT"],
                        "SER": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "OG", "HG", "C", "O", "OXT"],
                        "THR": ["N", "H", "CA", "HA", "CB", "HB", "CG2", "HG21", "HG22", "HG23", "OG1", "HG1", "C", "O", "OXT"],
                        "TRP": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "NE1", "HE1", "CE2", "CZ2", "HZ2", "CH2", "HH2", "CZ3", "HZ3", "CE3", "HE3", "CD2", "C", "O", "OXT"],
                        "TYR": ["N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "OH", "HH", "CE2", "HE2", "CD2", "HD2", "C", "O", "OXT"],
                        "VAL": ["N", "H", "CA", "HA", "CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23", "C", "O", "OXT"],
                        "NHE": ["N", "HN1", "HN2"],
                        "NME": ["N", "H", "CH3", "HH31", "HH32", "HH33"]
}
"""dict() that lists the atoms in the C-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the C-terminal version of that residue."""

RESIDUE_NTERMINAL_ATOM_LIST: Dict[str, List[str]] = {
                        'ACE': {'H1': ['CH3'], 'CH3': ['H1', 'H2', 'H3', 'C'], 'H2': ['CH3'], 'H3': ['CH3'], 'C': ['CH3', 'O'], 'O': ['C']},
                        'ALA': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB1', 'HB2', 'HB3'], 'HB1': ['CB'], 'HB2': ['CB'], 'HB3': ['CB'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'ARG': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG2', 'HG3', 'CD'], 'HG2': ['CG'], 'HG3': ['CG'], 'CD': ['CG', 'HD2', 'HD3', 'NE'], 'HD2': ['CD'], 'HD3': ['CD'], 'NE': ['CD', 'HE', 'CZ'], 'HE': ['NE'], 'CZ': ['NE', 'NH1', 'NH2'], 'NH1': ['CZ', 'HH11', 'HH12'], 'HH11': ['NH1'], 'HH12': ['NH1'], 'NH2': ['CZ', 'HH21', 'HH22'], 'HH21': ['NH2'], 'HH22': ['NH2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'ASN': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'OD1', 'ND2'], 'OD1': ['CG'], 'ND2': ['CG', 'HD21', 'HD22'], 'HD21': ['ND2'], 'HD22': ['ND2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'ASP': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'OD1', 'OD2'], 'OD1': ['CG'], 'OD2': ['CG'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'CYS': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'SG'], 'HB2': ['CB'], 'HB3': ['CB'], 'SG': ['CB', 'HG'], 'HG': ['SG'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'CYX': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'SG'], 'HB2': ['CB'], 'HB3': ['CB'], 'SG': ['CB'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'GLN': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG2', 'HG3', 'CD'], 'HG2': ['CG'], 'HG3': ['CG'], 'CD': ['CG', 'OE1', 'NE2'], 'OE1': ['CD'], 'NE2': ['CD', 'HE21', 'HE22'], 'HE21': ['NE2'], 'HE22': ['NE2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'GLU': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG2', 'HG3', 'CD'], 'HG2': ['CG'], 'HG3': ['CG'], 'CD': ['CG', 'OE1', 'OE2'], 'OE1': ['CD'], 'OE2': ['CD'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'GLY': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA2', 'HA3', 'C'], 'HA2': ['CA'], 'HA3': ['CA'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'HID': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'ND1', 'CD2'], 'ND1': ['CG', 'HD1', 'CE1'], 'HD1': ['ND1'], 'CE1': ['ND1', 'HE1', 'NE2'], 'HE1': ['CE1'], 'NE2': ['CE1', 'CD2'], 'CD2': ['CG', 'NE2', 'HD2'], 'HD2': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'HIE': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'ND1', 'CD2'], 'ND1': ['CG', 'CE1'], 'CE1': ['ND1', 'HE1', 'NE2'], 'HE1': ['CE1'], 'NE2': ['CE1', 'HE2', 'CD2'], 'HE2': ['NE2'], 'CD2': ['CG', 'NE2', 'HD2'], 'HD2': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'HIP': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'ND1', 'CD2'], 'ND1': ['CG', 'HD1', 'CE1'], 'HD1': ['ND1'], 'CE1': ['ND1', 'HE1', 'NE2'], 'HE1': ['CE1'], 'NE2': ['CE1', 'HE2', 'CD2'], 'HE2': ['NE2'], 'CD2': ['CG', 'NE2', 'HD2'], 'HD2': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'ILE': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB', 'CG2', 'CG1'], 'HB': ['CB'], 'CG2': ['CB', 'HG21', 'HG22', 'HG23'], 'HG21': ['CG2'], 'HG22': ['CG2'], 'HG23': ['CG2'], 'CG1': ['CB', 'HG12', 'HG13', 'CD1'], 'HG12': ['CG1'], 'HG13': ['CG1'], 'CD1': ['CG1', 'HD11', 'HD12', 'HD13'], 'HD11': ['CD1'], 'HD12': ['CD1'], 'HD13': ['CD1'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'LEU': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG', 'CD1', 'CD2'], 'HG': ['CG'], 'CD1': ['CG', 'HD11', 'HD12', 'HD13'], 'HD11': ['CD1'], 'HD12': ['CD1'], 'HD13': ['CD1'], 'CD2': ['CG', 'HD21', 'HD22', 'HD23'], 'HD21': ['CD2'], 'HD22': ['CD2'], 'HD23': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'LYS': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG2', 'HG3', 'CD'], 'HG2': ['CG'], 'HG3': ['CG'], 'CD': ['CG', 'HD2', 'HD3', 'CE'], 'HD2': ['CD'], 'HD3': ['CD'], 'CE': ['CD', 'HE2', 'HE3', 'NZ'], 'HE2': ['CE'], 'HE3': ['CE'], 'NZ': ['CE', 'HZ1', 'HZ2', 'HZ3'], 'HZ1': ['NZ'], 'HZ2': ['NZ'], 'HZ3': ['NZ'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'MET': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'HG2', 'HG3', 'SD'], 'HG2': ['CG'], 'HG3': ['CG'], 'SD': ['CG', 'CE'], 'CE': ['SD', 'HE1', 'HE2', 'HE3'], 'HE1': ['CE'], 'HE2': ['CE'], 'HE3': ['CE'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'PHE': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'CD1', 'CD2'], 'CD1': ['CG', 'HD1', 'CE1'], 'HD1': ['CD1'], 'CE1': ['CD1', 'HE1', 'CZ'], 'HE1': ['CE1'], 'CZ': ['CE1', 'HZ', 'CE2'], 'HZ': ['CZ'], 'CE2': ['CZ', 'HE2', 'CD2'], 'HE2': ['CE2'], 'CD2': ['CG', 'CE2', 'HD2'], 'HD2': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'PRO': {'N': ['H2', 'H3', 'CD', 'CA'], 'H2': ['N'], 'H3': ['N'], 'CD': ['N', 'HD2', 'HD3', 'CG'], 'HD2': ['CD'], 'HD3': ['CD'], 'CG': ['CD', 'HG2', 'HG3', 'CB'], 'HG2': ['CG'], 'HG3': ['CG'], 'CB': ['CG', 'HB2', 'HB3', 'CA'], 'HB2': ['CB'], 'HB3': ['CB'], 'CA': ['N', 'CB', 'HA', 'C'], 'HA': ['CA'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'SER': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'OG'], 'HB2': ['CB'], 'HB3': ['CB'], 'OG': ['CB', 'HG'], 'HG': ['OG'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'THR': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB', 'CG2', 'OG1'], 'HB': ['CB'], 'CG2': ['CB', 'HG21', 'HG22', 'HG23'], 'HG21': ['CG2'], 'HG22': ['CG2'], 'HG23': ['CG2'], 'OG1': ['CB', 'HG1'], 'HG1': ['OG1'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'TRP': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'CD1', 'CD2'], 'CD1': ['CG', 'HD1', 'NE1'], 'HD1': ['CD1'], 'NE1': ['CD1', 'HE1', 'CE2'], 'HE1': ['NE1'], 'CE2': ['NE1', 'CZ2', 'CD2'], 'CZ2': ['CE2', 'HZ2', 'CH2'], 'HZ2': ['CZ2'], 'CH2': ['CZ2', 'HH2', 'CZ3'], 'HH2': ['CH2'], 'CZ3': ['CH2', 'HZ3', 'CE3'], 'HZ3': ['CZ3'], 'CE3': ['CZ3', 'HE3', 'CD2'], 'HE3': ['CE3'], 'CD2': ['CG', 'CE2', 'CE3'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'TYR': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB2', 'HB3', 'CG'], 'HB2': ['CB'], 'HB3': ['CB'], 'CG': ['CB', 'CD1', 'CD2'], 'CD1': ['CG', 'HD1', 'CE1'], 'HD1': ['CD1'], 'CE1': ['CD1', 'HE1', 'CZ'], 'HE1': ['CE1'], 'CZ': ['CE1', 'OH', 'CE2'], 'OH': ['CZ', 'HH'], 'HH': ['OH'], 'CE2': ['CZ', 'HE2', 'CD2'], 'HE2': ['CE2'], 'CD2': ['CG', 'CE2', 'HD2'], 'HD2': ['CD2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']},
                        'VAL': {'N': ['H1', 'H2', 'H3', 'CA'], 'H1': ['N'], 'H2': ['N'], 'H3': ['N'], 'CA': ['N', 'HA', 'CB', 'C'], 'HA': ['CA'], 'CB': ['CA', 'HB', 'CG1', 'CG2'], 'HB': ['CB'], 'CG1': ['CB', 'HG11', 'HG12', 'HG13'], 'HG11': ['CG1'], 'HG12': ['CG1'], 'HG13': ['CG1'], 'CG2': ['CB', 'HG21', 'HG22', 'HG23'], 'HG21': ['CG2'], 'HG22': ['CG2'], 'HG23': ['CG2'], 'C': ['CA', 'O', '+1N'], 'O': ['C']}
}
"""dict() that lists the atoms in the N-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the N-terminal version of that residue."""

RESIDUE_CONNECTIVITY_MAP_NTERMINAL: Dict[str, Dict[str, List[str]]] = {
                        "ACE": {"H1": ["CH3"], "CH3": ["H1", "H2", "H3", "C"], "H2": ["CH3"], "H3": ["CH3"], "C": ["CH3", "O"], "O": ["C"]},
                        "ALA": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB1", "HB2", "HB3"], "HB1": ["CB"], "HB2": ["CB"], "HB3": ["CB"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "ARG": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "NE"], "HD2": ["CD"], "HD3": ["CD"], "NE": ["CD", "HE", "CZ"], "HE": ["NE"], "CZ": ["NE", "NH1", "NH2"], "NH1": ["CZ", "HH11", "HH12"], "HH11": ["NH1"], "HH12": ["NH1"], "NH2": ["CZ", "HH21", "HH22"], "HH21": ["NH2"], "HH22": ["NH2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "ASN": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "ND2"], "OD1": ["CG"], "ND2": ["CG", "HD21", "HD22"], "HD21": ["ND2"], "HD22": ["ND2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "ASP": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "OD1", "OD2"], "OD1": ["CG"], "OD2": ["CG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "CYS": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB", "HG"], "HG": ["SG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "CYX": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "SG"], "HB2": ["CB"], "HB3": ["CB"], "SG": ["CB"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "GLN": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "NE2"], "OE1": ["CD"], "NE2": ["CD", "HE21", "HE22"], "HE21": ["NE2"], "HE22": ["NE2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "GLU": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "OE1", "OE2"], "OE1": ["CD"], "OE2": ["CD"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "GLY": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA2", "HA3", "C"], "HA2": ["CA"], "HA3": ["CA"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "HID": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "CD2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "HIE": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "CE1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "HIP": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "ND1", "CD2"], "ND1": ["CG", "HD1", "CE1"], "HD1": ["ND1"], "CE1": ["ND1", "HE1", "NE2"], "HE1": ["CE1"], "NE2": ["CE1", "HE2", "CD2"], "HE2": ["NE2"], "CD2": ["CG", "NE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "ILE": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "CG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "CG1": ["CB", "HG12", "HG13", "CD1"], "HG12": ["CG1"], "HG13": ["CG1"], "CD1": ["CG1", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "LEU": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG", "CD1", "CD2"], "HG": ["CG"], "CD1": ["CG", "HD11", "HD12", "HD13"], "HD11": ["CD1"], "HD12": ["CD1"], "HD13": ["CD1"], "CD2": ["CG", "HD21", "HD22", "HD23"], "HD21": ["CD2"], "HD22": ["CD2"], "HD23": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "LYS": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "CD"], "HG2": ["CG"], "HG3": ["CG"], "CD": ["CG", "HD2", "HD3", "CE"], "HD2": ["CD"], "HD3": ["CD"], "CE": ["CD", "HE2", "HE3", "NZ"], "HE2": ["CE"], "HE3": ["CE"], "NZ": ["CE", "HZ1", "HZ2", "HZ3"], "HZ1": ["NZ"], "HZ2": ["NZ"], "HZ3": ["NZ"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "MET": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "HG2", "HG3", "SD"], "HG2": ["CG"], "HG3": ["CG"], "SD": ["CG", "CE"], "CE": ["SD", "HE1", "HE2", "HE3"], "HE1": ["CE"], "HE2": ["CE"], "HE3": ["CE"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "PHE": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "HZ", "CE2"], "HZ": ["CZ"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "PRO": {"N": ["H2", "H3", "CD", "CA"], "H2": ["N"], "H3": ["N"], "CD": ["N", "HD2", "HD3", "CG"], "HD2": ["CD"], "HD3": ["CD"], "CG": ["CD", "HG2", "HG3", "CB"], "HG2": ["CG"], "HG3": ["CG"], "CB": ["CG", "HB2", "HB3", "CA"], "HB2": ["CB"], "HB3": ["CB"], "CA": ["N", "CB", "HA", "C"], "HA": ["CA"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "SER": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "OG"], "HB2": ["CB"], "HB3": ["CB"], "OG": ["CB", "HG"], "HG": ["OG"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "THR": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG2", "OG1"], "HB": ["CB"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "OG1": ["CB", "HG1"], "HG1": ["OG1"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "TRP": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "NE1"], "HD1": ["CD1"], "NE1": ["CD1", "HE1", "CE2"], "HE1": ["NE1"], "CE2": ["NE1", "CZ2", "CD2"], "CZ2": ["CE2", "HZ2", "CH2"], "HZ2": ["CZ2"], "CH2": ["CZ2", "HH2", "CZ3"], "HH2": ["CH2"], "CZ3": ["CH2", "HZ3", "CE3"], "HZ3": ["CZ3"], "CE3": ["CZ3", "HE3", "CD2"], "HE3": ["CE3"], "CD2": ["CG", "CE2", "CE3"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "TYR": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB2", "HB3", "CG"], "HB2": ["CB"], "HB3": ["CB"], "CG": ["CB", "CD1", "CD2"], "CD1": ["CG", "HD1", "CE1"], "HD1": ["CD1"], "CE1": ["CD1", "HE1", "CZ"], "HE1": ["CE1"], "CZ": ["CE1", "OH", "CE2"], "OH": ["CZ", "HH"], "HH": ["OH"], "CE2": ["CZ", "HE2", "CD2"], "HE2": ["CE2"], "CD2": ["CG", "CE2", "HD2"], "HD2": ["CD2"], "C": ["CA", "O", "+1N"], "O": ["C"]},
                        "VAL": {"N": ["H1", "H2", "H3", "CA"], "H1": ["N"], "H2": ["N"], "H3": ["N"], "CA": ["N", "HA", "CB", "C"], "HA": ["CA"], "CB": ["CA", "HB", "CG1", "CG2"], "HB": ["CB"], "CG1": ["CB", "HG11", "HG12", "HG13"], "HG11": ["CG1"], "HG12": ["CG1"], "HG13": ["CG1"], "CG2": ["CB", "HG21", "HG22", "HG23"], "HG21": ["CG2"], "HG22": ["CG2"], "HG23": ["CG2"], "C": ["CA", "O", "+1N"], "O": ["C"]}
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the N-terminal version of each residue."""

RESIDUE_CATEGORIES: Dict[str, List[str]] = {
        "charged": ["R","H","K","D","E"],
        "negative": ["D","E"],
        "neutral": ["S","T","N","Q","C","Y","A","V","I","L","P","M","F","W","G"],
        "nonpolar": ["A","V","I","L","P","M","F","W","G"],
        "polar": ["R","H","K","D","E","S","T","N","Q","C","Y"],
        "positive": ["R","H","K"],
        "polar_neutral": ["S","T","N","Q","C","Y"],
        "strong_h_bond_donor": ["R","H","K","Y","C"],
        "weak_h_bond_donor": ["S","T","N","Q","W"],
}
"""dict() that describes basic characteristics of amino acids. Has (key,value) pairs
of ('characteric', list() of one-letter amino-acid codes). Covered characteristics are:
polar, nonpolar, charged, positive, negative, neutral, polar_neutral, strong_h_bond_donor,
weak_h_bond_donor."""

RESIDUE_VOLUME_MAPPER: Dict[str, float] = {
        "A": 88.6,
        "C": 108.5,
        "D": 111.1,
        "E": 138.4,
        "F": 189.9,
        "G": 60.1,
        "H": 153.2,
        "I": 166.7,
        "K": 168.6,
        "L": 166.7,
        "M": 162.9,
        "N": 114.1,
        "P": 112.7,
        "Q": 143.8,
        "R": 173.4,
        "S": 89.0,
        "T": 116.1,
        "V": 140.0,
        "W": 227.8,
        "Y": 193.6
}
"""dict() that maps one-letter amino-acid codes to their volume in cubic angstroms. 
source: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html"""

CAA_CHARGE_MAPPER: Dict[str, int] = {
        "ALA": 0,
        "CYS": 0,
        "ASP": -1,
        "GLU": -1,
        "PHE": 0,
        "GLY": 0,
        "ILE": 0,
        "LYS": 1,
        "LEU": 0,
        "MET": 0,
        "ASN": 0,
        "PRO": 0,
        "GLN": 0,
        "ARG": 1,
        "SER": 0,
        "THR": 0,
        "VAL": 0,
        "TRP": 0,
        "TYR": 0
}
"""dict() that maps three-letter canonical amino-acid codes to their formal charge. 
Note that HIS is temporarily removed from the map due to its potential 3 titration state.
This is used as the pool for target mutation selection when charge related keyword for mutation
is used"""

AA_CHARGE_MAPPER: Dict[str, int] = {
        "ALA": 0,
        "CYS": 0,
        "ASP": -1,
        "GLU": -1,
        "PHE": 0,
        "GLY": 0,
        "ILE": 0,
        "LYS": 1,
        "LEU": 0,
        "MET": 0,
        "ASN": 0,
        "PRO": 0,
        "GLN": 0,
        "ARG": 1,
        "SER": 0,
        "THR": 0,
        "VAL": 0,
        "TRP": 0,
        "TYR": 0,
        "HIS": 0,
        "HIE": 0,
        "HID": 0,
        "HIP": 1,
}
"""dict() that maps three-letter amino-acid codes to their formal charge. not limited to CAA.
This map includes residues of different titration states. This is used in determine the charge
of the original residue when charge related keyword for mutation is used"""

DEPROTONATION_MAPPER: Dict[str, Union[tuple, None]] = {
    "ASH": {"OD2": ("ASP", "HD2")},
    "CYS": {"SG": ("CYM", "HG")},
    "GLH": {"OE2": ("GLU", "HE2")},
    "HIP": {"ND1": ("HIE", "HD1"),"NE2": ("HID", "HE2")},
    "HID": {"ND1": ("HIE", "HD1")},  # note these are just switching
    "HIE": {"NE2": ("HID", "HE2")},  # note these are just switching
    "LYS": {"NZ": ("LYN", "HZ1")},
    "TYR": {"OH": ("TYM", "HH")},
    "ARG": {"NH2": ("AR0", "HH22"),"NH1": ("AR0", "HH12")}
}  # note this need to switch name after depro
"""a map of residue name of deprotonatable residues. This means the residue have ambiguous protonation state with a pH
around 7.0. The map is presented in the way that the deprotonatable residue is the key and its
target_atom : (deprotoned state, deprotoned atom) is the value. In the file, these residue name are arranged with each
charge state as an indent from +1 to -1.
{resi_name : {deprotonated_atom_name : (depro_resi_name, depro_proton_name), ...}
See /resource/ProtonationState.cdx for more detail"""

NOPROTON_LIST = ["ASP", "GLU", "MET"]
"""a list of residue name with no acidic proton"""

RESIDUE_NON_MUTATE_ATOM_MAPPER: Dict[str, List[str]] = {
    "default" : ["N", "H", "CA", "HA", "CB", "C", "O"],
    "PRO" : ['N','CA','HA','CB','C','O'],
    "GLY" : ['N','H','CA','C','O'],
}
"""a dictionary of atom names that won't change in a subsitution-mutation. Bascially mainchain+CB.
The key is the residue after the mutation."""

# yapf: enable

def convert_to_three_letter(one_letter: str) -> str:
    """Converts a one letter amino acid name to a three letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(one_letter) != 1:
        raise InvalidResidueCode(f"expecting one letter residue code. '{one_letter}' is invalid")
    one_letter = one_letter.upper()
    result = ONE_LETTER_AA_MAPPER.get(one_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {one_letter}")
    return result


def convert_to_canonical_three_letter(three_letter: str) -> str:
    """Converts a one letter amino acid name to a three letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(three_letter) != 3:
        raise InvalidResidueCode(f"expecting three letter residue code. '{three_letter}' is invalid")
    three_letter = three_letter.upper()
    result = THREE_TO_THREE_LETTER_CAA_MAPPER.get(three_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {three_letter}")
    return result


def convert_to_one_letter(three_letter: str) -> str:
    """Converts a threee letter amino acid name to a one letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(three_letter) != 3:
        raise InvalidResidueCode(f"expecting three letter residue code. '{three_letter}' is invalid")
    three_letter = three_letter.upper()
    result = THREE_LETTER_AA_MAPPER.get(three_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {three_letter}")
    return result


def get_element_aliases(ff: str, element: str) -> Set[str]:
    """Gets all element aliases for a given force field (ff) and element name, retungin in a set. If the ff is not supported, will log and exit."""
    if ff not in RESIDUE_ELEMENT_MAP:
        _LOGGER.error(f"{ff} is not a supported force field type. Allowed are '{', '.join(list(RESIDUE_ELEMENT_MAP.keys()))}'. Exiting...")
        exit(0)
    ff_dict = RESIDUE_ELEMENT_MAP[ff]
    aliases = []
    for alias, elem in ff_dict.items():
        if elem == element:
            aliases.append(alias)

    return set(aliases)


def one_letters_except(existing: str) -> List[str]:
    """Creates a list of all one letter amino acid codes. Case insensitive. If an invalid code is supplied, an enzy_htp.core.InvalidResidueCode() error is thrown."""
    existing = existing.upper()
    if existing not in ONE_LETTER_AA_MAPPER:
        raise InvalidResidueCode(f"{existing} is not a valid one letter amino acid")
    result = list(set(list(ONE_LETTER_AA_MAPPER.keys())))
    return list(filter(lambda s: s != existing, result))


def residue_polarity(code: str) -> str:
    """Determines the polarity of a one-letter nucleotide, being 'negative', 'neutral' or 'positive."""
    if len(code) != 1:
        raise InvalidResidueCode(f"expecting one letter residue code. '{code}' is invalid")

    result: str() = "unknown"
    for ptype in "positive negative neutral".split():
        if code in RESIDUE_CATEGORIES[ptype]:
            result = ptype
    return result


def non_polar(code: str) -> bool:
    # TODO(CJ): should probably check if it is a valid one letter residue code
    """Determines if a one-letter nucleotide amino acid is non-polar. Returns True if it is non-polar."""
    if len(code) != 1:
        raise InvalidResidueCode(f"expecting one letter residue code. '{code}' is invalid")

    return code in RESIDUE_CATEGORIES["nonpolar"]


def polar(code: str) -> bool:
    """Determines if a one-letter nucleotide amino acid is polar. Returns True if it is non-polar."""
    # TODO(CJ): should probably check if it is a valid one letter residue code
    return not non_polar(code)


def get_non_mutate_atom_names(residue_name: str) -> List[str]:
    """Get names of atoms that does not involve in a substitution mutation to
    the residue ({residue_name}). 
    For example: for mutating to most residues, all
    the side chain atoms will be replaced whiled the main chain atoms and CB won't
    change. So in this case, the returning list of atoms are the mainchain atoms and
    CB. For mutating to GLY, CB and HA are also removed so the list may change
    depending on the {residue_name}.
    This function is mainly used for mutate_stru_with_tleap().
    Args:
        residue_name: the 3-letter name for the residue after the mutation
    Returns:
        a list of atom names that remains the same before and after mutation"""
    if len(residue_name) != 3:
        raise InvalidResidueCode(f"expecting three letter residue code. '{residue_name}' is invalid")
    residue_name = residue_name.upper()
    result = RESIDUE_NON_MUTATE_ATOM_MAPPER.get(residue_name, RESIDUE_NON_MUTATE_ATOM_MAPPER["default"])

    return result

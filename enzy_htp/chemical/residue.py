"""Stores mappers and testing functions for residue information. Much of this functionality is later used in 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from typing import List, Set, Dict, Union
from ..core import InvalidResidueCode, _LOGGER


AA_LIST: List[str] = [
    "R",
    "H",
    "K",
    "D",
    "E",
    "S",
    "T",
    "N",
    "Q",
    "C",
    "G",
    "P",
    "A",
    "V",
    "I",
    "L",
    "M",
    "F",
    "Y",
    "W",
    "U"
]
"""Capitalized list of all one-letter amino acid names."""

THREE_LETTER_AA_MAPPER: Dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASH": "N",
    "ASN": "N",
    "ASP": "D",
    "CYM": "C",
    "CYS": "C",
    "CYX": "C",
    "GLH": "E",
    "GLN": "Q",
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
}
#TODO(shaoqz): @imp2 add related to canonical in the name
"""Contains mapping of all amino acids codes, with key value pairs of (three letter code, one letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_three_letter()"""

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
"""Contains mapping of all amino acids codes, with key value pairs of (one letter code, three letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_one_letter()"""

RESIDUE_ELEMENT_MAP: Dict[str, Dict[str, str]] = {
    "Amber": {
        "AG": "Ag",
        "AL": "Al",
        "BA": "Ba",
        "BE": "Be",
        "C": "C",
        "CA": "C",
        "CB": "C",
        "CD": "C",
        "CD1": "C",
        "CD2": "C",
        "CE": "C",
        "CE1": "C",
        "CE2": "C",
        "CE3": "C",
        "CG": "C",
        "CG1": "C",
        "CG2": "C",
        "CH2": "C",
        "CO": "Co",
        "CR": "Cr",
        "CS": "Cs",
        "CU": "Cu",
        "CZ": "C",
        "CZ2": "C",
        "CZ3": "C",
        "DY": "Dy",
        "ER": "Er",
        "EU": "Eu",
        "FE": "Fe",
        "GD": "Gd",
        "H": "H",
        "H1": "H",
        "H2": "H",
        "H3": "H",
        "HA": "H",
        "HA2": "H",
        "HA3": "H",
        "HB": "H",
        "HB1": "H",
        "HB2": "H",
        "HB3": "H",
        "HD1": "H",
        "HD11": "H",
        "HD12": "H",
        "HD13": "H",
        "HD2": "H",
        "HD21": "H",
        "HD22": "H",
        "HD23": "H",
        "HD3": "H",
        "HE": "H",
        "HE1": "H",
        "HE2": "H",
        "HE21": "H",
        "HE22": "H",
        "HE3": "H",
        "HF": "Hf",
        "HG": "H",
        "HG1": "H",
        "HG11": "H",
        "HG12": "H",
        "HG13": "H",
        "HG2": "H",
        "HG21": "H",
        "HG22": "H",
        "HG23": "H",
        "HG3": "H",
        "HH": "H",
        "HH11": "H",
        "HH12": "H",
        "HH2": "H",
        "HH21": "H",
        "HH22": "H",
        "HZ": "H",
        "HZ1": "H",
        "HZ2": "H",
        "HZ3": "H",
        "IN": "In",
        "K": "K",
        "LA": "La",
        "LI": "Li",
        "LU": "Lu",
        "MG": "Mg",
        "MN": "Mn",
        "N": "N",
        "NA": "Na",
        "ND": "Nd",
        "ND1": "N",
        "ND2": "N",
        "NE": "N",
        "NE1": "N",
        "NE2": "N",
        "NH1": "N",
        "NH2": "N",
        "NI": "Ni",
        "NZ": "N",
        "O": "O",
        "OD1": "O",
        "OD2": "O",
        "OE1": "O",
        "OE2": "O",
        "OG": "O",
        "OG1": "O",
        "OH": "O",
        "OXT": "O",
        "PB": "Pb",
        "PD": "Pd",
        "PR": "Pr",
        "PT": "Pt",
        "PU": "Pu",
        "RA": "Ra",
        "RB": "Rb",
        "SD": "S",
        "SG": "S",
        "SM": "Sm",
        "SN": "Sn",
        "SR": "Sr",
        "TB": "Tb",
        "TH": "Th",
        "TL": "Tl",
        "TM": "Tm",
        "U": "U",
        "V": "V",
        "Y": "Y",
        "YB": "Yb",
        "ZN": "Zn",
        "ZR": "Zr"
    }
}
"""Mapper that shows element names for alternative names of atoms. Key value structure is (force field name, (atom mapper)), where atom mapper is an additional
mapper that maps altered atom names to base atom names. Currently only defined for "Amber".

Example usage:

>>> initial_atom_name = "CA"
>>> RESIDUE_ELEMENT_MAP["Amber"][initial_atom_name] TODO also remove those ff keys
"C"

"""

RESIDUE_CONNECTIVITY_MAP: Dict[str, Dict[str, List[str]]] = {
    "ALA": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB1",
            "HB2",
            "HB3"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB1": [
            "CB"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "ARG": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "NE"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "CZ": [
            "NE",
            "NH1",
            "NH2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE": [
            "NE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HH11": [
            "NH1"
        ],
        "HH12": [
            "NH1"
        ],
        "HH21": [
            "NH2"
        ],
        "HH22": [
            "NH2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE": [
            "CD",
            "HE",
            "CZ"
        ],
        "NH1": [
            "CZ",
            "HH11",
            "HH12"
        ],
        "NH2": [
            "CZ",
            "HH21",
            "HH22"
        ],
        "O": [
            "C"
        ]
    },
    "ASH": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "OD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "OD2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ],
        "OD2": [
            "CG",
            "HD2"
        ]
    },
    "ASN": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "ND2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD21": [
            "ND2"
        ],
        "HD22": [
            "ND2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND2": [
            "CG",
            "HD21",
            "HD22"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ]
    },
    "ASP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "OD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ],
        "OD2": [
            "CG"
        ]
    },
    "CYM": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB3",
            "HB2",
            "SG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "SG": [
            "CB"
        ]
    },
    "CYS": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "SG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "SG": [
            "CB",
            "HG"
        ]
    },
    "CYX": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "SG": [
            "CB"
        ]
    },
    "GLH": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "OE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE2": [
            "OE2"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ],
        "OE2": [
            "CD",
            "HE2"
        ]
    },
    "GLN": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE21": [
            "NE2"
        ],
        "HE22": [
            "NE2"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE2": [
            "CD",
            "HE21",
            "HE22"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ]
    },
    "GLU": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "OE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ],
        "OE2": [
            "CD"
        ]
    },
    "GLY": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA2",
            "HA3",
            "C"
        ],
        "H": [
            "N"
        ],
        "HA2": [
            "CA"
        ],
        "HA3": [
            "CA"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "HID": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "HIE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "HIP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "HYP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "CB",
            "HA",
            "C"
        ],
        "CB": [
            "CG",
            "HB2",
            "HB3",
            "CA"
        ],
        "CD": [
            "N",
            "HD22",
            "HD23",
            "CG"
        ],
        "CG": [
            "CD",
            "HG",
            "OD1",
            "CB"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "OD1"
        ],
        "HD22": [
            "CD"
        ],
        "HD23": [
            "CD"
        ],
        "HG": [
            "CG"
        ],
        "N": [
            "CD",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG",
            "HD1"
        ]
    },
    "ILE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "CG1"
        ],
        "CD1": [
            "CG1",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CG1": [
            "CB",
            "HG12",
            "HG13",
            "CD1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "LEU": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CD2": [
            "CG",
            "HD21",
            "HD22",
            "HD23"
        ],
        "CG": [
            "CB",
            "HG",
            "CD1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HD21": [
            "CD2"
        ],
        "HD22": [
            "CD2"
        ],
        "HD23": [
            "CD2"
        ],
        "HG": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "LYN": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "CE"
        ],
        "CE": [
            "CD",
            "HE2",
            "HE3",
            "NZ"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HZ2": [
            "NZ"
        ],
        "HZ3": [
            "NZ"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NZ": [
            "CE",
            "HZ2",
            "HZ3"
        ],
        "O": [
            "C"
        ]
    },
    "LYS": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "CE"
        ],
        "CE": [
            "CD",
            "HE2",
            "HE3",
            "NZ"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HZ1": [
            "NZ"
        ],
        "HZ2": [
            "NZ"
        ],
        "HZ3": [
            "NZ"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NZ": [
            "CE",
            "HZ1",
            "HZ2",
            "HZ3"
        ],
        "O": [
            "C"
        ]
    },
    "MET": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CE": [
            "SD",
            "HE1",
            "HE2",
            "HE3"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "SD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE1": [
            "CE"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "SD": [
            "CG",
            "CE"
        ]
    },
    "PHE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "HZ",
            "CE2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HZ": [
            "CZ"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "PRO": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "CB",
            "HA",
            "C"
        ],
        "CB": [
            "CG",
            "HB2",
            "HB3",
            "CA"
        ],
        "CD": [
            "N",
            "HD2",
            "HD3",
            "CG"
        ],
        "CG": [
            "CD",
            "HG2",
            "HG3",
            "CB"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "CD",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "SER": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "OG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "OG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OG": [
            "CB",
            "HG"
        ]
    },
    "THR": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "OG1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG1": [
            "OG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OG1": [
            "CB",
            "HG1"
        ]
    },
    "TRP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "NE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "CE3"
        ],
        "CE2": [
            "NE1",
            "CZ2",
            "CD2"
        ],
        "CE3": [
            "CZ3",
            "HE3",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CH2": [
            "CZ2",
            "HH2",
            "CZ3"
        ],
        "CZ2": [
            "CE2",
            "HZ2",
            "CH2"
        ],
        "CZ3": [
            "CH2",
            "HZ3",
            "CE3"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HE1": [
            "NE1"
        ],
        "HE3": [
            "CE3"
        ],
        "HH2": [
            "CH2"
        ],
        "HZ2": [
            "CZ2"
        ],
        "HZ3": [
            "CZ3"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE1": [
            "CD1",
            "HE1",
            "CE2"
        ],
        "O": [
            "C"
        ]
    },
    "TYR": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "OH",
            "CE2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HH": [
            "OH"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OH": [
            "CZ",
            "HH"
        ]
    },
    "VAL": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG1",
            "CG2"
        ],
        "CG1": [
            "CB",
            "HG11",
            "HG12",
            "HG13"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG11": [
            "CG1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ]
    },
    "WAT": {
        "H1": [
            "O"
        ],
        "H2": [
            "O"
        ],
        "O": [
            "H1",
            "H2"
        ]
    }
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms)."""

RESIDUE_CONNECTIVITY_MAP_CTERMINAL: Dict[str, Dict[str, List[str]]] = {
    "ALA": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB1",
            "HB2",
            "HB3"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB1": [
            "CB"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "ARG": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "NE"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "CZ": [
            "NE",
            "NH1",
            "NH2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE": [
            "NE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HH11": [
            "NH1"
        ],
        "HH12": [
            "NH1"
        ],
        "HH21": [
            "NH2"
        ],
        "HH22": [
            "NH2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE": [
            "CD",
            "HE",
            "CZ"
        ],
        "NH1": [
            "CZ",
            "HH11",
            "HH12"
        ],
        "NH2": [
            "CZ",
            "HH21",
            "HH22"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "ASN": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "ND2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD21": [
            "ND2"
        ],
        "HD22": [
            "ND2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND2": [
            "CG",
            "HD21",
            "HD22"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ],
        "OXT": [
            "C"
        ]
    },
    "ASP": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "OD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ],
        "OD2": [
            "CG"
        ],
        "OXT": [
            "C"
        ]
    },
    "CYS": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "SG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ],
        "SG": [
            "CB",
            "HG"
        ]
    },
    "CYX": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ],
        "SG": [
            "CB"
        ]
    },
    "GLN": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE21": [
            "NE2"
        ],
        "HE22": [
            "NE2"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE2": [
            "CD",
            "HE21",
            "HE22"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ],
        "OXT": [
            "C"
        ]
    },
    "GLU": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "OE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ],
        "OE2": [
            "CD"
        ],
        "OXT": [
            "C"
        ]
    },
    "GLY": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA2",
            "HA3",
            "C"
        ],
        "H": [
            "N"
        ],
        "HA2": [
            "CA"
        ],
        "HA3": [
            "CA"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "HID": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "CD2"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "HIE": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "HIP": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "HYP": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "CB",
            "HA",
            "C"
        ],
        "CB": [
            "CG",
            "HB2",
            "HB3",
            "CA"
        ],
        "CD": [
            "N",
            "HD22",
            "HD23",
            "CG"
        ],
        "CG": [
            "CD",
            "HG",
            "OD1",
            "CB"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "OD1"
        ],
        "HD22": [
            "CD"
        ],
        "HD23": [
            "CD"
        ],
        "HG": [
            "CG"
        ],
        "N": [
            "CD",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG",
            "HD1"
        ],
        "OXT": [
            "C"
        ]
    },
    "ILE": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "CG1"
        ],
        "CD1": [
            "CG1",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CG1": [
            "CB",
            "HG12",
            "HG13",
            "CD1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "LEU": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CD2": [
            "CG",
            "HD21",
            "HD22",
            "HD23"
        ],
        "CG": [
            "CB",
            "HG",
            "CD1",
            "CD2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HD21": [
            "CD2"
        ],
        "HD22": [
            "CD2"
        ],
        "HD23": [
            "CD2"
        ],
        "HG": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "LYS": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "CE"
        ],
        "CE": [
            "CD",
            "HE2",
            "HE3",
            "NZ"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HZ1": [
            "NZ"
        ],
        "HZ2": [
            "NZ"
        ],
        "HZ3": [
            "NZ"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NZ": [
            "CE",
            "HZ1",
            "HZ2",
            "HZ3"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "MET": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CE": [
            "SD",
            "HE1",
            "HE2",
            "HE3"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "SD"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE1": [
            "CE"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ],
        "SD": [
            "CG",
            "CE"
        ]
    },
    "NHE": {
        "HN1": [
            "N"
        ],
        "HN2": [
            "N"
        ],
        "N": [
            "HN1",
            "HN2",
            "-1C"
        ]
    },
    "NME": {
        "CH3": [
            "N",
            "HH31",
            "HH32",
            "HH33"
        ],
        "H": [
            "N"
        ],
        "HH31": [
            "CH3"
        ],
        "HH32": [
            "CH3"
        ],
        "HH33": [
            "CH3"
        ],
        "N": [
            "H",
            "CH3",
            "-1C"
        ]
    },
    "PHE": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "HZ",
            "CE2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HZ": [
            "CZ"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "PRO": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "CB",
            "HA",
            "C"
        ],
        "CB": [
            "CG",
            "HB2",
            "HB3",
            "CA"
        ],
        "CD": [
            "N",
            "HD2",
            "HD3",
            "CG"
        ],
        "CG": [
            "CD",
            "HG2",
            "HG3",
            "CB"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "CD",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "SER": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "OG"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "OG"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OG": [
            "CB",
            "HG"
        ],
        "OXT": [
            "C"
        ]
    },
    "THR": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "OG1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG1": [
            "OG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OG1": [
            "CB",
            "HG1"
        ],
        "OXT": [
            "C"
        ]
    },
    "TRP": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "NE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "CE3"
        ],
        "CE2": [
            "NE1",
            "CZ2",
            "CD2"
        ],
        "CE3": [
            "CZ3",
            "HE3",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CH2": [
            "CZ2",
            "HH2",
            "CZ3"
        ],
        "CZ2": [
            "CE2",
            "HZ2",
            "CH2"
        ],
        "CZ3": [
            "CH2",
            "HZ3",
            "CE3"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HE1": [
            "NE1"
        ],
        "HE3": [
            "CE3"
        ],
        "HH2": [
            "CH2"
        ],
        "HZ2": [
            "CZ2"
        ],
        "HZ3": [
            "CZ3"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "NE1": [
            "CD1",
            "HE1",
            "CE2"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    },
    "TYR": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "OH",
            "CE2"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HH": [
            "OH"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OH": [
            "CZ",
            "HH"
        ],
        "OXT": [
            "C"
        ]
    },
    "VAL": {
        "C": [
            "CA",
            "O",
            "OXT"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG1",
            "CG2"
        ],
        "CG1": [
            "CB",
            "HG11",
            "HG12",
            "HG13"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG11": [
            "CG1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H",
            "CA",
            "-1C"
        ],
        "O": [
            "C"
        ],
        "OXT": [
            "C"
        ]
    }
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the C-terminal version of each residue."""

RESIDUE_CTERMINAL_ATOM_LIST: Dict[str, List[str]] = {
    "ACE": [
        "H1",
        "CH3",
        "H2",
        "H3",
        "C",
        "O"
    ],
    "ALA": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB1",
        "HB2",
        "HB3",
        "C",
        "O"
    ],
    "ARG": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG2",
        "HG3",
        "CD",
        "HD2",
        "HD3",
        "NE",
        "HE",
        "CZ",
        "NH1",
        "HH11",
        "HH12",
        "NH2",
        "HH21",
        "HH22",
        "C",
        "O"
    ],
    "ASN": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "OD1",
        "ND2",
        "HD21",
        "HD22",
        "C",
        "O"
    ],
    "ASP": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "OD1",
        "OD2",
        "C",
        "O"
    ],
    "CYS": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "SG",
        "HG",
        "C",
        "O"
    ],
    "CYX": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "SG",
        "C",
        "O"
    ],
    "GLN": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG2",
        "HG3",
        "CD",
        "OE1",
        "NE2",
        "HE21",
        "HE22",
        "C",
        "O"
    ],
    "GLU": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG2",
        "HG3",
        "CD",
        "OE1",
        "OE2",
        "C",
        "O"
    ],
    "GLY": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA2",
        "HA3",
        "C",
        "O"
    ],
    "HID": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "ND1",
        "HD1",
        "CE1",
        "HE1",
        "NE2",
        "CD2",
        "HD2",
        "C",
        "O"
    ],
    "HIE": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "ND1",
        "CE1",
        "HE1",
        "NE2",
        "HE2",
        "CD2",
        "HD2",
        "C",
        "O"
    ],
    "HIP": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "ND1",
        "HD1",
        "CE1",
        "HE1",
        "NE2",
        "HE2",
        "CD2",
        "HD2",
        "C",
        "O"
    ],
    "ILE": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB",
        "CG2",
        "HG21",
        "HG22",
        "HG23",
        "CG1",
        "HG12",
        "HG13",
        "CD1",
        "HD11",
        "HD12",
        "HD13",
        "C",
        "O"
    ],
    "LEU": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG",
        "CD1",
        "HD11",
        "HD12",
        "HD13",
        "CD2",
        "HD21",
        "HD22",
        "HD23",
        "C",
        "O"
    ],
    "LYS": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG2",
        "HG3",
        "CD",
        "HD2",
        "HD3",
        "CE",
        "HE2",
        "HE3",
        "NZ",
        "HZ1",
        "HZ2",
        "HZ3",
        "C",
        "O"
    ],
    "MET": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "HG2",
        "HG3",
        "SD",
        "CE",
        "HE1",
        "HE2",
        "HE3",
        "C",
        "O"
    ],
    "PHE": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "CD1",
        "HD1",
        "CE1",
        "HE1",
        "CZ",
        "HZ",
        "CE2",
        "HE2",
        "CD2",
        "HD2",
        "C",
        "O"
    ],
    "PRO": [
        "N",
        "H2",
        "H3",
        "CD",
        "HD2",
        "HD3",
        "CG",
        "HG2",
        "HG3",
        "CB",
        "HB2",
        "HB3",
        "CA",
        "HA",
        "C",
        "O"
    ],
    "SER": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "OG",
        "HG",
        "C",
        "O"
    ],
    "THR": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB",
        "CG2",
        "HG21",
        "HG22",
        "HG23",
        "OG1",
        "HG1",
        "C",
        "O"
    ],
    "TRP": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "CD1",
        "HD1",
        "NE1",
        "HE1",
        "CE2",
        "CZ2",
        "HZ2",
        "CH2",
        "HH2",
        "CZ3",
        "HZ3",
        "CE3",
        "HE3",
        "CD2",
        "C",
        "O"
    ],
    "TYR": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB2",
        "HB3",
        "CG",
        "CD1",
        "HD1",
        "CE1",
        "HE1",
        "CZ",
        "OH",
        "HH",
        "CE2",
        "HE2",
        "CD2",
        "HD2",
        "C",
        "O"
    ],
    "VAL": [
        "N",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA",
        "CB",
        "HB",
        "CG1",
        "HG11",
        "HG12",
        "HG13",
        "CG2",
        "HG21",
        "HG22",
        "HG23",
        "C",
        "O"
    ]
}
"""dict() that lists the atoms in the C-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the C-terminal version of that residue."""

RESIDUE_NTERMINAL_ATOM_LIST: Dict[str, List[str]] = {
   'ACE': ['H1', 'CH3', 'H2', 'H3', 'C', 'O'],
                        'ALA': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O'],
                        'ARG': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O'],
                        'ASN': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O'],
                        'ASP': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2', 'C', 'O'],
                        'CYS': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'HG', 'C', 'O'],
                        'CYX': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG', 'C', 'O'],
                        'GLN': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O'],
                        'GLU': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'C', 'O'],
                        'GLY': ['N', 'H1', 'H2', 'H3', 'CA', 'HA2', 'HA3', 'C', 'O'],
                        'HID': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2', 'C', 'O'],
                        'HIE': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
                        'HIP': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
                        'ILE': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'C', 'O'],
                        'LEU': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O'],
                        'LYS': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'],
                        'MET': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O'],
                        'PHE': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
                        'PRO': ['N', 'H2', 'H3', 'CD', 'HD2', 'HD3', 'CG', 'HG2', 'HG3', 'CB', 'HB2', 'HB3', 'CA', 'HA', 'C', 'O'],
                        'SER': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'OG', 'HG', 'C', 'O'],
                        'THR': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'OG1', 'HG1', 'C', 'O'],
                        'TRP': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O'],
                        'TYR': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
                        'VAL': ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O'] }
"""dict() that lists the atoms in the N-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the N-terminal version of that residue."""

RESIDUE_CONNECTIVITY_MAP_NTERMINAL: Dict[str, Dict[str, List[str]]] = {
    "ACE": {
        "C": [
            "CH3",
            "O"
        ],
        "CH3": [
            "H1",
            "H2",
            "H3",
            "C"
        ],
        "H1": [
            "CH3"
        ],
        "H2": [
            "CH3"
        ],
        "H3": [
            "CH3"
        ],
        "O": [
            "C"
        ]
    },
    "ALA": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB1",
            "HB2",
            "HB3"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB1": [
            "CB"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "ARG": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "NE"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "CZ": [
            "NE",
            "NH1",
            "NH2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE": [
            "NE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HH11": [
            "NH1"
        ],
        "HH12": [
            "NH1"
        ],
        "HH21": [
            "NH2"
        ],
        "HH22": [
            "NH2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "NE": [
            "CD",
            "HE",
            "CZ"
        ],
        "NH1": [
            "CZ",
            "HH11",
            "HH12"
        ],
        "NH2": [
            "CZ",
            "HH21",
            "HH22"
        ],
        "O": [
            "C"
        ]
    },
    "ASN": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "ND2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD21": [
            "ND2"
        ],
        "HD22": [
            "ND2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "ND2": [
            "CG",
            "HD21",
            "HD22"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ]
    },
    "ASP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CG": [
            "CB",
            "OD1",
            "OD2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "OD1": [
            "CG"
        ],
        "OD2": [
            "CG"
        ]
    },
    "CYS": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "SG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "SG": [
            "CB",
            "HG"
        ]
    },
    "CYX": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "SG"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "SG": [
            "CB"
        ]
    },
    "GLN": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE21": [
            "NE2"
        ],
        "HE22": [
            "NE2"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "NE2": [
            "CD",
            "HE21",
            "HE22"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ]
    },
    "GLU": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "OE1",
            "OE2"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "OE1": [
            "CD"
        ],
        "OE2": [
            "CD"
        ]
    },
    "GLY": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA2",
            "HA3",
            "C"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA2": [
            "CA"
        ],
        "HA3": [
            "CA"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "HID": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "HIE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "ND1": [
            "CG",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "HIP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD2": [
            "CG",
            "NE2",
            "HD2"
        ],
        "CE1": [
            "ND1",
            "HE1",
            "NE2"
        ],
        "CG": [
            "CB",
            "ND1",
            "CD2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "ND1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "NE2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "ND1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "NE2": [
            "CE1",
            "HE2",
            "CD2"
        ],
        "O": [
            "C"
        ]
    },
    "ILE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "CG1"
        ],
        "CD1": [
            "CG1",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CG1": [
            "CB",
            "HG12",
            "HG13",
            "CD1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "LEU": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD11",
            "HD12",
            "HD13"
        ],
        "CD2": [
            "CG",
            "HD21",
            "HD22",
            "HD23"
        ],
        "CG": [
            "CB",
            "HG",
            "CD1",
            "CD2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD11": [
            "CD1"
        ],
        "HD12": [
            "CD1"
        ],
        "HD13": [
            "CD1"
        ],
        "HD21": [
            "CD2"
        ],
        "HD22": [
            "CD2"
        ],
        "HD23": [
            "CD2"
        ],
        "HG": [
            "CG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "LYS": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD": [
            "CG",
            "HD2",
            "HD3",
            "CE"
        ],
        "CE": [
            "CD",
            "HE2",
            "HE3",
            "NZ"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "CD"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "HZ1": [
            "NZ"
        ],
        "HZ2": [
            "NZ"
        ],
        "HZ3": [
            "NZ"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "NZ": [
            "CE",
            "HZ1",
            "HZ2",
            "HZ3"
        ],
        "O": [
            "C"
        ]
    },
    "MET": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CE": [
            "SD",
            "HE1",
            "HE2",
            "HE3"
        ],
        "CG": [
            "CB",
            "HG2",
            "HG3",
            "SD"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HE1": [
            "CE"
        ],
        "HE2": [
            "CE"
        ],
        "HE3": [
            "CE"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "SD": [
            "CG",
            "CE"
        ]
    },
    "PHE": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "HZ",
            "CE2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HZ": [
            "CZ"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "PRO": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "CB",
            "HA",
            "C"
        ],
        "CB": [
            "CG",
            "HB2",
            "HB3",
            "CA"
        ],
        "CD": [
            "N",
            "HD2",
            "HD3",
            "CG"
        ],
        "CG": [
            "CD",
            "HG2",
            "HG3",
            "CB"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD2": [
            "CD"
        ],
        "HD3": [
            "CD"
        ],
        "HG2": [
            "CG"
        ],
        "HG3": [
            "CG"
        ],
        "N": [
            "H2",
            "H3",
            "CD",
            "CA"
        ],
        "O": [
            "C"
        ]
    },
    "SER": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "OG"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HG": [
            "OG"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "OG": [
            "CB",
            "HG"
        ]
    },
    "THR": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG2",
            "OG1"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG1": [
            "OG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "OG1": [
            "CB",
            "HG1"
        ]
    },
    "TRP": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "NE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "CE3"
        ],
        "CE2": [
            "NE1",
            "CZ2",
            "CD2"
        ],
        "CE3": [
            "CZ3",
            "HE3",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CH2": [
            "CZ2",
            "HH2",
            "CZ3"
        ],
        "CZ2": [
            "CE2",
            "HZ2",
            "CH2"
        ],
        "CZ3": [
            "CH2",
            "HZ3",
            "CE3"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HE1": [
            "NE1"
        ],
        "HE3": [
            "CE3"
        ],
        "HH2": [
            "CH2"
        ],
        "HZ2": [
            "CZ2"
        ],
        "HZ3": [
            "CZ3"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "NE1": [
            "CD1",
            "HE1",
            "CE2"
        ],
        "O": [
            "C"
        ]
    },
    "TYR": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB2",
            "HB3",
            "CG"
        ],
        "CD1": [
            "CG",
            "HD1",
            "CE1"
        ],
        "CD2": [
            "CG",
            "CE2",
            "HD2"
        ],
        "CE1": [
            "CD1",
            "HE1",
            "CZ"
        ],
        "CE2": [
            "CZ",
            "HE2",
            "CD2"
        ],
        "CG": [
            "CB",
            "CD1",
            "CD2"
        ],
        "CZ": [
            "CE1",
            "OH",
            "CE2"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB2": [
            "CB"
        ],
        "HB3": [
            "CB"
        ],
        "HD1": [
            "CD1"
        ],
        "HD2": [
            "CD2"
        ],
        "HE1": [
            "CE1"
        ],
        "HE2": [
            "CE2"
        ],
        "HH": [
            "OH"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ],
        "OH": [
            "CZ",
            "HH"
        ]
    },
    "VAL": {
        "C": [
            "CA",
            "O",
            "+1N"
        ],
        "CA": [
            "N",
            "HA",
            "CB",
            "C"
        ],
        "CB": [
            "CA",
            "HB",
            "CG1",
            "CG2"
        ],
        "CG1": [
            "CB",
            "HG11",
            "HG12",
            "HG13"
        ],
        "CG2": [
            "CB",
            "HG21",
            "HG22",
            "HG23"
        ],
        "H1": [
            "N"
        ],
        "H2": [
            "N"
        ],
        "H3": [
            "N"
        ],
        "HA": [
            "CA"
        ],
        "HB": [
            "CB"
        ],
        "HG11": [
            "CG1"
        ],
        "HG12": [
            "CG1"
        ],
        "HG13": [
            "CG1"
        ],
        "HG21": [
            "CG2"
        ],
        "HG22": [
            "CG2"
        ],
        "HG23": [
            "CG2"
        ],
        "N": [
            "H1",
            "H2",
            "H3",
            "CA"
        ],
        "O": [
            "C"
        ]
    }
}
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the N-terminal version of each residue."""


RESIDUE_CATEGORIES: Dict[str, List[str]] = {
    "charged": [
        "R",
        "H",
        "K",
        "D",
        "E"
    ],
    "negative": [
        "D",
        "E"
    ],
    "neutral": [
        "S",
        "T",
        "N",
        "Q",
        "C",
        "Y",
        "A",
        "V",
        "I",
        "L",
        "P",
        "M",
        "F",
        "W",
        "G"
    ],
    "nonpolar": [
        "A",
        "V",
        "I",
        "L",
        "P",
        "M",
        "F",
        "W",
        "G"
    ],
    "polar": [
        "R",
        "H",
        "K",
        "D",
        "E",
        "S",
        "T",
        "N",
        "Q",
        "C",
        "Y"
    ],
    "positive": [
        "R",
        "H",
        "K"
    ]
}
"""dict() that describes basic characteristics of amino acids. Has (key,value) pairs
of ('characteric', list() of one-letter amino-acid codes). Covered characteristics are:
polar, nonpolar, charged, positive, negative, neutral."""

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
source: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html
"""

DEPROTONATION_MAPPER: Dict[str, Union[tuple, None]] = { "ASH" : {"OD2" : ("ASP","HD2")},
                                                        "CYS" : {"SG" : ("CYM","HG" )},
                                                        "GLH" : {"OE2" : ("GLU","HE2")},
                                                "HIP" :{"ND1" : ("HIE", "HD1"),
                                                        "NE2" : ("HID", "HE2")},
                                                "HID" : {"ND1" : ("HIE", "HD1")}, # note these are just switching
                                                "HIE" : {"NE2" : ("HID", "HE2")}, # note these are just switching
                                                "LYS" : {"NZ" : ("LYN","HZ1")},
                                                        "TYR" : {"OH" : ("TYM","HH")},
                                                "ARG" : {"NH2" : ("AR0","HH22"),
                                                         "NH1" : ("AR0","HH12")}} # note this need to switch name after depro
"""a map of residue name of deprotonatable residues. This means the residue have ambiguous protonation state with a pH
around 7.0. The map is presented in the way that the deprotonatable residue is the key and its
target_atom : (deprotoned state, deprotoned atom) is the value. In the file, these residue name are arranged with each
charge state as an indent from +1 to -1.
{resi_name : {deprotonated_atom_name : (depro_resi_name, depro_proton_name), ...}
See /resource/ProtonationState.cdx for more detail"""

NOPROTON_LIST = ["ASP", "GLU", "MET"]
"""a list of residue name with no acidic proton"""


def convert_to_three_letter(one_letter: str) -> str:
    """Converts a one letter amino acid name to a three letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(one_letter) != 1:
        raise InvalidResidueCode(
            f"expecting one letter residue code. '{one_letter}' is invalid"
        )
    one_letter = one_letter.upper()
    result = ONE_LETTER_AA_MAPPER.get(one_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {one_letter}")
    return result


def convert_to_one_letter(three_letter: str) -> str:
    """Converts a threee letter amino acid name to a one letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(three_letter) != 3:
        raise InvalidResidueCode(
            f"expecting three letter residue code. '{three_letter}' is invalid"
        )
    three_letter = three_letter.upper()
    result = THREE_LETTER_AA_MAPPER.get(three_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {three_letter}")
    return result


def get_element_aliases(ff: str, element: str) -> Set[str]:
    """Gets all element aliases for a given force field (ff) and element name, retungin in a set. If the ff is not supported, will log and exit."""
    if ff not in RESIDUE_ELEMENT_MAP:
        _LOGGER.error(
            f"{ff} is not a supported force field type. Allowed are '{', '.join(list(RESIDUE_ELEMENT_MAP.keys()))}'. Exiting..."
        )
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
        raise InvalidResidueCode(
            f"expecting one letter residue code. '{code}' is invalid"
        )

    result: str() = "unknown"
    for ptype in "positive negative neutral".split():
        if code in RESIDUE_CATEGORIES[ptype]:
            result = ptype
    return result


def non_polar(code: str) -> bool:
    # TODO(CJ): should probably check if it is a valid one letter residue code
    """Determines if a one-letter nucleotide amino acid is non-polar. Returns True if it is non-polar."""
    if len(code) != 1:
        raise InvalidResidueCode(
            f"expecting one letter residue code. '{code}' is invalid"
        )

    return code in RESIDUE_CATEGORIES["nonpolar"]


def polar(code: str) -> bool:
    """Determines if a one-letter nucleotide amino acid is polar. Returns True if it is non-polar."""
    # TODO(CJ): should probably check if it is a valid one letter residue code
    return not non_polar(code)

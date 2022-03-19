"""Stores mappers and testing functions for residue information. Much of this functionality is later used in 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
"""

from typing import List, Set, Dict

from ..core import InvalidResidueCode, _LOGGER


AA_LIST : List[str]= [
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
    "U",
]
"""Capitalized list of all one-letter amino acid names."""

THREE_LETTER_AA_MAPPER : Dict[str,str] = {
    "ARG": "R",
    "HIS": "H",
    "HIE": "H",
    "HIP": "H",
    "HID": "H",
    "LYS": "K",
    "ASP": "D",
    "GLU": "E",
    "SER": "S",
    "THR": "T",
    "ASN": "N",
    "GLN": "Q",
    "CYS": "C",
    "SEC": "U",
    "GLY": "G",
    "PRO": "P",
    "ALA": "A",
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "MET": "M",
    "PHE": "F",
    "TYR": "Y",
    "TRP": "W",
}
"""Contains mapping of all amino acids codes, with key value pairs of (three letter code, one letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_three_letter()"""

ONE_LETTER_AA_MAPPER : Dict[str,str] = {
    "R": "ARG",
    "H": "HIS",
    "K": "LYS",
    "D": "ASP",
    "E": "GLU",
    "S": "SER",
    "T": "THR",
    "N": "ASN",
    "Q": "GLN",
    "C": "CYS",
    "U": "SEC",
    "G": "GLY",
    "P": "PRO",
    "A": "ALA",
    "V": "VAL",
    "I": "ILE",
    "L": "LEU",
    "M": "MET",
    "F": "PHE",
    "Y": "TYR",
    "W": "TRP",
}
"""Contains mapping of all amino acids codes, with key value pairs of (one letter code, three letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_one_letter()"""

RESIDUE_ELEMENT_MAP : Dict[str,Dict[str,str]] = {
    "Amber": {
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
        "CZ": "C",
        "CZ2": "C",
        "CZ3": "C",
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
        "N": "N",
        "ND1": "N",
        "ND2": "N",
        "NE": "N",
        "NE1": "N",
        "NE2": "N",
        "NH1": "N",
        "NH2": "N",
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
        "SD": "S",
        "SG": "S",
        "LI": "Li",
        "NA": "Na",
        "K": "K",
        "RB": "Rb",
        "CS": "Cs",
        "MG": "Mg",
        "TL": "Tl",
        "CU": "Cu",
        "AG": "Ag",
        "BE": "Be",
        "NI": "Ni",
        "PT": "Pt",
        "ZN": "Zn",
        "CO": "Co",
        "PD": "Pd",
        "CR": "Cr",
        "FE": "Fe",
        "V": "V",
        "MN": "Mn",
        "YB": "Yb",
        "SN": "Sn",
        "PB": "Pb",
        "EU": "Eu",
        "SR": "Sr",
        "SM": "Sm",
        "BA": "Ba",
        "RA": "Ra",
        "AL": "Al",
        "IN": "In",
        "Y": "Y",
        "LA": "La",
        "PR": "Pr",
        "ND": "Nd",
        "GD": "Gd",
        "TB": "Tb",
        "DY": "Dy",
        "ER": "Er",
        "TM": "Tm",
        "LU": "Lu",
        "HF": "Hf",
        "ZR": "Zr",
        "U": "U",
        "PU": "Pu",
        "TH": "Th",
    }
}
"""Mapper that shows element names for alternative names of atoms. Key value structure is (force field name, (atom mapper)), where atom mapper is an additional
mapper that maps altered atom names to base atom names. Currently only defined for "Amber".

Example usage:

>>> initial_atom_name = "CA"
>>> RESIDUE_ELEMENT_MAP["Amber"][initial_atom_name]
"C"

"""

def convert_to_three_letter(one_letter: str) -> str:
    """Converts a one letter amino acid name to a three letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(one_letter) != 1:
        raise InvalidResidueCode(
            f"expecting one letter residue code. '{one_letter}' is invalid")
    one_letter = one_letter.upper()
    result = ONE_LETTER_AA_MAPPER.get(one_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {one_letter}")
    return result


def convert_to_one_letter(three_letter: str) -> str:
    """Converts a threee letter amino acid name to a one letter. If supplied code is invalid, raises an enzy_htp.core.InvalidResidueCode() exception."""
    if len(three_letter) != 3:
        raise InvalidResidueCode(
            f"expecting three letter residue code. '{three_letter}' is invalid")
    three_letter = three_letter.upper()
    result = THREE_LETTER_AA_MAPPER.get(three_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {three_letter}")
    return result


def get_element_aliases(ff: str, element: str) -> Set[str]:
    """Gets all element aliases for a given force field (ff) and element name, retungin in a set. If the ff is not supported, will log and exit."""
    if ff not in ff_dict:
        _LOGGER.error(f"{ff} is not a supported force field type. Allowed are '{', '.join(list(RESIDUE_ELEMENT_MAP.keys()))}'. Exiting...")
        exit( 0 )
    ff_dict = RESIDUE_ELEMENT_MAP[ff]
    aliases = []
    for alias, elem in ff_dict.items():
        if elem == element:
            aliases.append(alias)

    return set(aliases)


def one_letters_except(existing: str) -> List[str]:
    """Creates a list of all one letter amino acid codes. If an invalid code is supplied, an enzy_htp.core.InvalidResidueCode() error is thrown."""
    if existing not in ONE_LETTER_AA_MAPPER:
        raise InvalidResidueCode(f"{existing} is not a valid one letter amino acid")
    result = list(ONE_LETTER_AA_MAPPER.keys())
    return list(filter(lambda s: s != existing, result))

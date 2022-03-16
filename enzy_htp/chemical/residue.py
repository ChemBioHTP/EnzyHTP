"""TODO: DOCUMENATION"""

from ..core import InvalidResidueCode

# TODO maybe check for non-canonical?

THREE_LETTER_AA_MAPPER = {
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


ONE_LETTER_AA_MAPPER = {
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


def convert_to_three_letter(one_letter: str) -> str:
    """TODO DOCUMENTATION"""
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
    """TODO DOCUMENTATION"""
    if len(three_letter) != 3:
        raise InvalidResidueCode(
            f"expecting three letter residue code. '{three_letter}' is invalid"
        )
    three_letter = three_letter.upper()
    result = THREE_LETTER_AA_MAPPER.get(three_letter, None)
    if not result:
        raise InvalidResidueCode(f"Invalid residue code {three_letter}")
    return result

"""Stores mappers and definitions for different types of metals often found in PDBs.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
"""

import sys
from typing import Dict, Union

from enzy_htp.core import _LOGGER

METAL_MAPPER: Dict[str, str] = {
    "LI": "Li",
    "NA": "Na",
    "Na+": "Na",
    "K": "K",
    "K+": "K",
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
    "HG": "Hg",
    "CD": "Cd",
    "YB": "Yb",
    "CA": "Ca",
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
    "CE": "Ce",
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

METAL_CENTER_MAP: Dict[str, str] = {
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
    "HG": "Hg",
    "CD": "Cd",
    "YB": "Yb",
    "CA": "Ca",
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
    "CE": "Ce",
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
"""Contains mapping of metal types. Used to check if a residue atom name from a PDB is a metal."""

IONIC_RADII: Dict[str, Union[int, None]] = {
    "N": 1.32,
    "O": 1.26,
    "S": 1.70,
    "Mg": 0.86,
    "Li": 0.90,
    "Zn": 0.88,
    "Co": 0.68,  # +3 ls
    "Fe": 0.92,  # +2 hs / PDB disctance based (5HIO) / can not handel really well
    "Mn": 0.90,  # +2 / modified based on PDB (1K20)
    "Ca": 1.14,
    "Cu": 0.87,  # +2
    "Hg": 1.10,
}
"""Mapping of metal elements to ionic radii. Value is 'None' if non-existent. Reference: doi:10.1107/S0567739476001551 """

VDW_RADII = {
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "Mg": 1.73,
    "Li": 1.82,
    "Tl": 1.96,
    "Cu": 1.40,
    "Ag": 1.72,
    "Be": 1.53,
    "Ni": 1.63,
    "Pt": 1.75,
    "Zn": 1.39,
    "Pd": 1.63,
    "Hg": 1.55,
    "Cd": 1.58,
    "Ca": 2.31,
    "Sn": 2.17,
    "Pb": 2.02,
    "Sr": 2.49,
    "Ba": 2.68,
    "Ra": 2.83,
    "Al": 1.84,
    "In": 1.93,
}
"""Mapping of metal elements to Van-Der Waals (VDW) radii. Value is 'None' if non-existent."""

DONOR_ATOM_LIST = [
    "NH2", "NE", "NH1", "ND1", "NE2", "NZ", "OD1", 'OD2', "OE1", 'OE2', "OG", "OG1", "ND2", "OD1", "OE1", "NE2", "SG", "SD", "OH", "NE1"
]
"""List for atom names that can be a qaulify electron donor atom to a coordination center
The dictionary key here was for the parsing logic of atom names. Current the PDB format
TODO see if there are more format and the way to decouple this from atom.name. (maybe just
use the most popular one parse the rest into it. In that way there's no need for these dict
keys)
See /resource/AtomName.cdx for more detail"""


def get_atom_radii(element: str, method: str = "ionic") -> float:
    """Given a element name and a method, returns the radii for that element. Allowed methods are 'ionic' and 'vdw'.
    Exits if an invalid radii type is given or no record found of this element in the mapper"""
    if method == "ionic":
        result = IONIC_RADII.get(element, None)
    elif method == "vdw":
        result = VDW_RADII.get(element, None)
    else:
        _LOGGER.error(f"The radii method '{method}' is not allowed. Only 'ionic' and 'vdw' are supported. Exiting...")
        sys.exit(1)
    if result is None:
        _LOGGER.error(f"query element {element} not in {method} mapper. Consider add it in.")
        sys.exit(1)
    return result


# TODO(CJ): add method is_metal() that checks if an atom is a metal

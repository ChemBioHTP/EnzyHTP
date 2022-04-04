"""Stores mappers and definitions for different types of metals often found in PDBs.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
"""

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


IONIC_RADII : Dict[str,Union[int,None]] =  {
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
}
"""Mapping of metal elements to ionic radii. Value is 'None' if non-existent. Reference: doi:10.1107/S0567739476001551 """

VDW_RADII =  {
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

def get_metal_radii( atom_name : str, method : str = 'ionic' ) -> Union[float,None]:
    """Given a metal atom and a method, returns the radii for that metal type. Allowed methods are 'ionic' and 'vdw'. Exits if an invalid radii type is given. 
	Returns 'None' if not defined."""
    #TODO(CJ): provide some basic input atom sanitization. (i.e. LI -> Li, li -> Li)
    if method == 'ionic':
        return IONIC_RADII.get( atom_name, None )
    elif method == 'vdw':
        return VDW_RADII.get( atom_name, None )
    else:
        _LOGGER.error(f"The radii method '{method}' is not allowed. Only 'ionic' and 'vdw' are supported. Exiting...")
        exit( 1 )
     

#TODO(CJ): add method is_metal() that checks if an atom is a metal

"""
chemical module for enzy_htp. Responsibilities include identification and mapping of residues, metals and solvents
as well as providing infrastructure to store relevanta chemical information in a databse at enzy_htp/chemical/chemical-database.json.


Author: Chris Jurich <chris.jurich@vanerbilt.edu>
Date: 2022-03-19
"""
from .atoms import get_h_bond_length
from .enum import ResidueType
from .metal import METAL_MAPPER, METAL_CENTER_MAP, get_metal_radii
from .solvent import RD_SOLVENT_LIST, RD_NON_LIGAND_LIST
from .residue import (
    THREE_LETTER_AA_MAPPER,
    ONE_LETTER_AA_MAPPER,
    RESIDUE_VOLUME_MAPPER,
    RESIDUE_CATEGORIES,
    one_letters_except,
    convert_to_one_letter,
    convert_to_three_letter,
    get_element_aliases,
    residue_polarity,
    non_polar,
)

"""
knowledge module for enzy_htp. This is all the basic science knowledge database 
of EnzyHTP. Responsibilities include identification and mapping of residues, metals and solvents


Author: Chris Jurich <chris.jurich@vanerbilt.edu>
Date: 2022-03-19
"""
from .atoms import get_h_bond_length, get_valid_generic_atom_name
from .enum import ResidueType
from .metal import METAL_MAPPER, METAL_CENTER_MAP, get_atom_radii
from .solvent import RD_SOLVENT_LIST, RD_NON_LIGAND_LIST
from .residue import (
    THREE_LETTER_AA_MAPPER,
    ONE_LETTER_AA_MAPPER,
    RESIDUE_VOLUME_MAPPER,
    RESIDUE_CATEGORIES,
    AA_CHARGE_MAPPER,
    one_letters_except,
    convert_to_one_letter,
    convert_to_three_letter,
    get_element_aliases,
    residue_polarity,
    non_polar,
)

from .enzyme_commission import parse_ec_number

from .level_of_theory import (
    LevelOfTheory,
    QMLevelOfTheory,
    MMLevelOfTheory,
)

from .physics import electric_field_strength

"""Stores mappers and testing functions for residue information. Much of this functionality is later used in 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from typing import List, Set, Dict
from .db import load_from_db
from ..core import InvalidResidueCode, _LOGGER


AA_LIST: List[str] = load_from_db("AA_LIST")
"""Capitalized list of all one-letter amino acid names."""

THREE_LETTER_AA_MAPPER: Dict[str, str] = load_from_db("THREE_LETTER_AA_MAPPER") #@shaoqz: @imp2 add related to canonical in the name
"""Contains mapping of all amino acids codes, with key value pairs of (three letter code, one letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_three_letter()"""

ONE_LETTER_AA_MAPPER: Dict[str, str] = load_from_db("ONE_LETTER_AA_MAPPER")
"""Contains mapping of all amino acids codes, with key value pairs of (one letter code, three letter code). Should NOT be called directly for code conversion. Instead used enzy_htp.chemical.residue.convert_to_one_letter()"""

RESIDUE_ELEMENT_MAP: Dict[str, Dict[str, str]] = load_from_db("RESIDUE_ELEMENT_MAP")
"""Mapper that shows element names for alternative names of atoms. Key value structure is (force field name, (atom mapper)), where atom mapper is an additional
mapper that maps altered atom names to base atom names. Currently only defined for "Amber".

Example usage:

>>> initial_atom_name = "CA"
>>> RESIDUE_ELEMENT_MAP["Amber"][initial_atom_name]
"C"

"""

RESIDUE_CONNECTIVITY_MAP: Dict[str, Dict[str, List[str]]] = load_from_db(
    "RESIDUE_CONNECTIVITY_MAP"
)
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms)."""

RESIDUE_CONNECTIVITY_MAP_CTERMINAL: Dict[str, Dict[str, List[str]]] = load_from_db(
    "RESIDUE_CONNECTIVITY_MAP_CTERMINAL"
)
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the C-terminal version of each residue."""

RESIDUE_CTERMINAL_ATOM_LIST: Dict[str, List[str]] = load_from_db(
    "RESIDUE_CTERMINAL_ATOM_LIST"
)
"""dict() that lists the atoms in the C-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the C-terminal version of that residue."""

RESIDUE_CTERMINAL_ATOM_LIST: Dict[str, List[str]] = load_from_db(
    "RESIDUE_CTERMINAL_ATOM_LIST"
)
"""dict() that lists the atoms in the N-terminal version of a residue. (key, value) pairs are
(3-letter AA name, list of atom names) where the list of atom names are the modified list of 
atom names for the N-terminal version of that residue."""

RESIDUE_CONNECTIVITY_MAP_NTERMINAL: Dict[str, Dict[str, List[str]]] = load_from_db(
    "RESIDUE_CONNECTIVITY_MAP_NTERMINAL"
)
"""dict() that maps the connectivity of atom's based on their parent residue and atom identity.
(key, value) pairs are (3-letter AA name, connector), where connector is another dict() with 
(key, value) pairs of (atom name, list of connected atoms). Similar to RESIDUE_CONNECTIVITY_MAP
except has the mappings for the N-terminal version of each residue."""


RESIDUE_CATEGORIES: Dict[str, List[str]] = load_from_db("RESIDUE_CATEGORIES")
"""dict() that describes basic characteristics of amino acids. Has (key,value) pairs
of ('characteric', list() of one-letter amino-acid codes). Covered characteristics are:
polar, nonpolar, charged, positive, negative, neutral."""

RESIDUE_VOLUME_MAPPER: Dict[str, float] = load_from_db("RESIDUE_VOLUME_MAPPER")
"""dict() that maps one-letter amino-acid codes to their volume in cubic angstroms. 
source: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html
"""


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

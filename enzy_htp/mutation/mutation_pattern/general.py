"""Functions related to decode the mutation_pattern. The pattern indicates a set of mutations.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""
from typing import List
import re

from enzy_htp.core.exception import InvalidMutationPatternSyntax
from enzy_htp.structure import Structure
from ..mutation import Mutation


def decode_mutation_pattern(stru: Structure, pattern: str) -> List[Mutation]:
    """
    decode the mutation {pattern} and return a list of mutation_obj using {stru} as the reference

    The mutation pattern is defined below:
    *Pattern Syntax:*
        "section_1,section_2,section_3,..."
        The pattern is composed by sections seperate by comma and each section can be one of the
        format below:
        1. direct indication            : XA###Y
        2. random N mutation in a set   : r:N[mutation_set_patterns]
        3. all mutation in a set: a     : a:[mutation_set_patterns]

        The mutation_set_patterns is seperated by comma and each describes 2 things:
        1. position_pattern: a set of positions
                            (using the selection syntax in the selection module)
        2. target_aa_pattern: a set of target mutations apply to all positions in the current set
                            (using syntax in the target_aa_pattern module)
        The two pattern are seperated by ":" and a mutation_set_patterns looks like:
        "position_pattern_0:target_aa_pattern_0, ..."

        Overall an example of pattern will be:
        "RA154W, r:100[resi 289 around 4 and not resi 36:larger, proj(id 1000, id 2023, positive, 10):more_negative_charge]"
        * here proj() is a hypothetical selection function
    """
    sections = seperate_sections(pattern)
    for section in sections:
        section_type = get_section_type(section)

def seperate_sections(pattern:str) -> List[str]:
    """seperate a mutation pattern into sections
    *Note that current method use re and do not support [] in a []
    (e.g. 1,2,3,a:[4,[5,6],7] does not work)"""

    section_seperate_pattern = r"[^,[]+(?:\[[^\]]*\])*[^,]*"
    sections = re.findall(section_seperate_pattern, pattern.strip())
    sections = [i.strip() for i in sections]

    return sections

def get_section_type(section_pattern: str) -> str:
    """determine the type of the mutation section base on the pattern
    the supported types now are:
    + d: direct indication
    + r: random over the set
    + a: all mutations in the set"""
    if section_pattern.startswith("r"):
        return "r"
    if section_pattern.startswith("a"):
        return "a"
    if section_pattern[0].isupper():
        return "d"
    raise InvalidMutationPatternSyntax(
        f"Mutation pattern section startswith unsupported letter {section_pattern[0]}. Cant determine type."
        )

def decode_direct_mutation():
    """decode the mutation pattern section that directly indicate the mutation.
    Return a list of mutation objects."""

def decode_random_mutation():
    """decode the mutation pattern section that random over the mutation set.
    Return a list of mutation objects."""

def decode_all_mutation():
    """decode the mutation pattern section that mutate all in the mutation set.
    Return a list of mutation objects."""


TYPE_SECTION_DECODERS = {
    "d" : decode_direct_mutation,
    "r" : decode_random_mutation,
    "a" : decode_all_mutation,
}

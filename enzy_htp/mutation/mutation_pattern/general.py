"""Functions related to decode the mutation_pattern. The pattern indicates a set of mutations.
Main API: decode_mutation_pattern

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""
from typing import Dict, List
import re

from enzy_htp.core.exception import InvalidMutationPatternSyntax
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Structure
from ..mutation import Mutation, decode_mutation_flag, is_valid_mutation


def decode_mutation_pattern(stru: Structure, pattern: str) -> List[Mutation]:
    """
    decode the mutation {pattern} and return a list of mutation_obj using {stru} as the reference
    The syntax of the mutation pattern is defined in `enzy_htp.mutation.general.assign_mutation`
    """
    result_mutants = []
    mutant_patterns = seperate_mutant_patterns(pattern)
    for mutant_pattern in mutant_patterns:
        mutant_pattern = mutant_pattern.strip("{}")
        p_mutant_mapper = {}
        section_patterns = seperate_section_patterns(mutant_pattern)
        for section_pattern in section_patterns:
            section_type = get_section_type(section_pattern)
            # here we decode it into a list of mutants each composed by a list of Mutation
            p_mutant_mapper[section_pattern] = TYPE_SECTION_DECODERS[section_type](stru, section_pattern)
        pattern_mutants = combine_section_mutation(p_mutant_mapper)
        result_mutants.extend(pattern_mutants)
    return result_mutants

def seperate_mutant_patterns(pattern:str) -> List[str]:
    """seperate a mutation pattern into pattern of each mutants"""
    section_seperate_pattern = r"(?:[^,[{]|(?:\{[^\}]*\})|(?:\[[^\]]*\]))+[^,]*"
    mutants = re.findall(section_seperate_pattern, pattern.strip())
    mutants = [i.strip() for i in mutants]

    return mutants

def seperate_section_patterns(pattern:str) -> List[str]:
    """seperate a mutation pattern into sections
    *Note that current method use re and do not support [] in a []
    (e.g. 1,2,3,a:[4,[5,6],7] does not work)"""

    section_seperate_pattern = r"(?:[^,[{]|(?:\[[^\]]*\]))+[^,]*"
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

def decode_direct_mutation(stru: Structure, section_pattern: str) -> List[Mutation]:
    """decode the mutation pattern section that directly indicate the mutation.
    Return a list of mutation objects."""
    mutation_obj = decode_mutation_flag(section_pattern)
    is_valid_mutation(mutation_obj, stru)
    return mutation_obj

def decode_random_mutation(stru: Structure, section_pattern: str) -> List[Mutation]:
    """decode the mutation pattern section that random over the mutation set.
    Return a list of mutation objects."""
    result = ""
    return result

def decode_all_mutation(stru: Structure, section_pattern: str) -> List[Mutation]:
    """decode the mutation pattern section that mutate all in the mutation set.
    Return a list of mutation objects."""
    result = ""
    return result


TYPE_SECTION_DECODERS = {
    "d" : decode_direct_mutation,
    "r" : decode_random_mutation,
    "a" : decode_all_mutation,
}

def combine_section_mutation(mutation_mapper: Dict[str, Mutation]) -> List[List[Mutation]]:
    """Combine mutations decoded from each section. Return a list of final mutants"""
"""Functions related to decode the mutation_pattern. The pattern indicates a set of mutations.
Main API: decode_mutation_pattern

Author: QZ Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-01-26
"""
from typing import Dict, List, Set, Tuple
import re
import itertools

from enzy_htp.core.exception import InvalidMutationPatternSyntax
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import (get_random_list_elem, pop_random_list_elem, product_lists_allow_empty, split_but_brackets)
from enzy_htp.structure import Structure
from enzy_htp.mutation_class import (
    Mutation, 
    generate_from_mutation_flag, 
    generate_mutation_from_target_list,
)
from .position_pattern import decode_position_pattern
from .target_aa_pattern import check_target_aa_pattern, decode_target_aa_pattern


def decode_mutation_pattern(stru: Structure, pattern: str) -> List[List[Mutation]]:
    """
    decode the mutation {pattern} and return a list of mutants (list of mutation_obj)
    using {stru} as the reference. The syntax of the mutation pattern is defined in
    `enzy_htp.mutation.api.assign_mutation`
    Args:
        stru: the reference structure
        pattern: the mutation pattern to be decode
    Return:
        A list of mutants, each is a list of Mutation()
    """
    result_mutants = []
    mutant_patterns = seperate_mutant_patterns(pattern)
    for mutant_pattern in mutant_patterns:
        mutant_pattern = mutant_pattern.strip("{}")
        p_mutant_mapper: Dict[str, List[List[Mutation]]] = {}
        section_patterns = seperate_section_patterns(mutant_pattern)
        for section_pattern in section_patterns:
            section_type = get_section_type(section_pattern)
            # here we decode it into a list of mutants each composed by a list of Mutation
            p_mutant_mapper[section_pattern] = TYPE_SECTION_DECODERS[section_type](stru, section_pattern)
        pattern_mutants = combine_section_mutant(p_mutant_mapper)
        result_mutants.extend(pattern_mutants)
    return result_mutants


def seperate_mutant_patterns(pattern: str) -> List[str]:
    """seperate a mutation pattern into pattern of each mutants"""
    section_seperate_pattern = r"(?:[^,[{]|(?:\{[^\}]*\})|(?:\[[^\]]*\]))+[^,]*" #TODO use the one from general
    mutants = re.findall(section_seperate_pattern, pattern.strip())
    mutants = [i.strip() for i in mutants]

    return mutants


def seperate_section_patterns(pattern: str) -> List[str]:
    """seperate a mutation pattern into sections
    *Note that current method use re and do not support [] in a []
    (e.g. 1,2,3,a:[4,[5,6],7] does not work)"""

    section_seperate_pattern = r"(?:[^,[{]|(?:\[[^\]]*\]))+[^,]*" #TODO use the one from general
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
    raise InvalidMutationPatternSyntax(f"Mutation pattern section startswith unsupported letter {section_pattern[0]}. Cant determine type.")


# d:
def decode_direct_mutation(stru: Structure, section_pattern: str) -> List[List[Mutation]]:
    """decode the mutation pattern section that directly indicate the mutation.
    Return a list of mutation objects.
    pattern_example: XA###Y"""
    mutation_obj = generate_from_mutation_flag(section_pattern)
    mutation_obj.is_valid_mutation(stru)
    return [[mutation_obj]]


# r:
def decode_random_mutation(stru: Structure, section_pattern: str) -> List[List[Mutation]]:
    """decode the mutation pattern section that random over the mutation set.
    Return a list of mutation objects. (M number of N point mutants)
    Args:
        stru: the reference structure
        section_pattern: the section pattern to be decode
    Return:
        A list of mutants, each is a list of Mutation()

    pattern_example:
        r:N[xxx:yyy]*M or r:NR[xxx:yyy]*MR
        (M number of N point mutants)
        (R stands for whether repeating mutation is allowed for each mutant and
        whether repeating mutant is allowed in the result, respectively)"""
    re_pattern = r"r:([0-9]*)(R?)\[(.+)\]\*([0-9]*)(R?)"
    mut_point_num, point_allow_repeat, mutation_esm_patterns, mutant_num, mutant_allow_repeat = re.match(re_pattern,
                                                                                                         section_pattern).groups()
    mut_point_num = int(mut_point_num)
    mutant_num = int(mutant_num)
    mutation_esm_mapper = decode_mutation_esm_pattern(stru, mutation_esm_patterns)  # {mutation_site: Mutation}

    _LOGGER.info(f"generating random mutants in positions: {list(mutation_esm_mapper.keys())} ({len(mutation_esm_mapper)} sites total)")
    if len(mutation_esm_mapper) < mut_point_num:
        raise InvalidMutationPatternSyntax(
            f"number of desired point mutations are more than the total number of possible mutation sites in the ensemble, desired: {mut_point_num}, possible_sites: {len(mutation_esm_mapper)}"
        )

    result: List[List[Mutation]] = []
    while len(result) < mutant_num:
        # I. make mutations for each mutant
        each_mutant: Dict[tuple, Mutation] = {}  # point mutation of each mutant

        if not point_allow_repeat:
            non_repeat_points = list(mutation_esm_mapper.keys())

        temp_point_num = mut_point_num  # will only not eq to mut_point_num when repeat is allowed
        while len(each_mutant) < temp_point_num:
            # 1. determine positionn
            if point_allow_repeat:
                new_position = get_random_list_elem(list(mutation_esm_mapper.keys()))
                if new_position in each_mutant:
                    _LOGGER.warning(
                        f"repeating mutation is generated for {new_position}, the later one is used, also less num of mutations in this mutant (point_allow_repeat: True)"
                    )
                    temp_point_num -= 1  # control the number of mutations
            else:  # point_allow_repeat is None
                new_position = pop_random_list_elem(non_repeat_points)
            # 2. determine target
            new_mutation = get_random_list_elem(mutation_esm_mapper[new_position])
            each_mutant[new_position] = new_mutation  # use dict to make sure each point only have 1 mutation #yapf: disable

        # II. consider if repeating or not and add valid mutant
        if mutant_allow_repeat:
            result.append(set(each_mutant.values()))
        else:  # do not allow repeat mutant -- check if repeat
            each_mutant = set(each_mutant.values())
            if each_mutant in result:
                _LOGGER.info(f"repeating MUTANT is generated: {each_mutant}, regenerating it. (point_allow_repeat: False)")
                continue
            result.append(each_mutant)

    result = [list(x) for x in result]
    return result


# a:
def decode_all_mutation(stru: Structure, section_pattern: str) -> List[List[Mutation]]:
    """decode the mutation pattern section that mutate all in the mutation set.
    There will be the maxium same number of mutations as the number of positions.
    If "M" is not specificed, non-mutation of each site will also be included.
    e.g.: 'a_site' : ILM, 'b_site' : ILM will give you 16 mutants containing 2 point,
    1 point and a WT.
    If "M" is specificed, each point will be mutated and there will be the same number of
    point mutations in each mutant as the number of positions in the mutation_esm_pattern.
    If "X" is specificed, only <=X points will be mutated. 
    (e.g.: X=2 will only generate double, single point mutations, and WT)
    If both "X" and "M" are specified, strictly X points will be mutated. 
    (e.g.: X=1 will only generate single point mutations.)
    Args:
        stru: the reference structure
        section_pattern: the section pattern to be decode
    Return:
        A list of mutants, each is a list of Mutation()

    pattern_example:
        a:[xxx:yyy] or a:M[xxx:yyy] or a:X[xxx:yyy]"""
    re_pattern = r"a:(\d*)(M?)\[(.+)\]"
    mutation_cap_num, force_mutate_each_point, mutation_esm_patterns = re.match(
        re_pattern, section_pattern
        ).groups()
    mutation_esm_mapper = decode_mutation_esm_pattern(stru, mutation_esm_patterns)  #{position: mutations}
    # mutation cap: X
    target_positions: List[List] = [list(mutation_esm_mapper.keys()),]
    if mutation_cap_num:
        mutation_cap_num = int(mutation_cap_num)
        if mutation_cap_num < len(mutation_esm_mapper):
            # This case we need new target positions as combinations of picking X from all positions
            target_positions = itertools.combinations(target_positions[0], mutation_cap_num)
    
    result = set()
    for positions in target_positions:
        position_mutation_esm = [mutation_esm_mapper[i] for i in positions]
        # the force mutate flag: M
        if force_mutate_each_point:
            posi_result = list(itertools.product(*position_mutation_esm))
        else:
            posi_result = product_lists_allow_empty(position_mutation_esm)
        for mut in posi_result:
            result.add(tuple(mut))
    
    # remove duplicated ones 
    # (since when X < total_sites and not force mutate there will be repeating ones)
    result = [list(i) for i in result]

    return result


TYPE_SECTION_DECODERS = {
    "d": decode_direct_mutation,
    "r": decode_random_mutation,
    "a": decode_all_mutation,
}


def decode_mutation_esm_pattern(stru: Structure, mutation_esm_patterns: str) -> Dict[Tuple[str, int], List[Mutation]]:
    """decode mutation esm pattern into a list of mutation objects
    pattern_example: position_pattern_0:target_aa_pattern_0, ...
    NOTE: if different mutation_esm_pattern have shared points.
          The target of the shared points with be the union of all
    Args:
        stru: the reference structure
        mutation_esm_patterns: the mutation_esm_pattern to be decode
    Return:
        A dictionary of position : mutant (a list of Mutation()) pair

    pattern_example:
        position:target_aa"""

    esm_result: Dict[tuple, List[Mutation]] = {}
    esm_pattern_list = [i.strip() for i in split_but_brackets(mutation_esm_patterns, ",")]
    for esm_pattern in esm_pattern_list:
        position_pattern, target_aa_pattern = esm_pattern.split(":")
        check_target_aa_pattern(target_aa_pattern)  # give warning about pattern
        esm_positions = decode_position_pattern(stru, position_pattern, if_name=True)
        for esm_position, orig_resi in esm_positions:
            posi_target_aa = decode_target_aa_pattern(orig_resi, target_aa_pattern)
            posi_mutation = generate_mutation_from_target_list(esm_position, orig_resi, posi_target_aa)
            if esm_position in esm_result:  # shared position case
                esm_result[esm_position] = list(set(esm_result[esm_position]) | set(posi_mutation))
            else:
                esm_result[esm_position] = posi_mutation
    # logging
    _LOGGER.info("Mutation ensemble of the section:")
    for k, v in esm_result.items():
        _LOGGER.info(f"    {k} : {[x.get_target(if_one_letter=True) for x in v]} (total: {len(v)})")
    _LOGGER.info(f"(total: {len(esm_result)} position)")

    return esm_result


def combine_section_mutant(mutation_mapper: Dict[str, List[List[Mutation]]]) -> List[List[Mutation]]:
    """Combine mutants decoded from each section. Return a list of final mutants.
    Args:
        mutation_mapper: {section_str : list of mutant}
    NOTE: merge repeating mutations into one and result permutation of combination
          of each section"""
    if len(mutation_mapper) == 1:
        return list(mutation_mapper.values())[0]  # save a lot of time
    return [list(set(itertools.chain.from_iterable(x))) for x in itertools.product(*mutation_mapper.values())]

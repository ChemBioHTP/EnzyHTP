"""Submodule that defines the Mutation namedtuple() which describes a single point mutation in an 
enzyme. Additionally provides utility functions for determnining if Mutation()'s are defined
and if they satisfy certain change requirements.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-15
"""
import re
import string
from collections import namedtuple
from typing import List, Dict, Tuple

from enzy_htp.core.exception import InvalidMutationFlagSyntax
from enzy_htp.core.logger import _LOGGER
from enzy_htp.chemical import ONE_LETTER_AA_MAPPER
import enzy_htp.chemical as chem
import enzy_htp.structure as es

Mutation = namedtuple("Mutation", "orig target chain_id res_idx")
Mutation.__doc__ = f"""Named tuple representing a single point mutation in an enzyme.
   

	Attributes:
   		orig: the one-letter code of the original amino acid. Can be in [ "{", ".join(ONE_LETTER_AA_MAPPER.keys())}"].
		target: the one-letter code of the target mutation. Can be in [ "{", ".join(ONE_LETTER_AA_MAPPER.keys())}"].
		chain_id: a single capital letter.
		res_idx: the 1-indexed int() of the residue to Mutate

        *In the case of WT, the tuple is defined as (None, "WT", None, None)
"""

def decode_mutation_flag(mutation_flag: str) -> Mutation:
    """XA##Y -> ("X", "Y", "A", ##)
    WT -> (None, "WT", None, None)
    XA##X -> (None, "WT", None, None)
    *we may need to support 3-letter mutation in the future like: TYQA##TYP"""

    mutation_flag = mutation_flag.strip()
    if mutation_flag == "WT":
        return Mutation(None, "WT", None, None)
    pattern = r'([A-Z])([A-Z])?([0-9]+)([A-Z])'
    flag_match = re.match(pattern, mutation_flag)
    if flag_match is None:
        raise InvalidMutationFlagSyntax(
            f"{mutation_flag} doesnt match ([A-Z])([A-Z])?([0-9]+)([A-Z])"
            )

    orig = flag_match.group(1)
    chain_id = flag_match.group(2)
    res_idx = flag_match.group(3)
    target = flag_match.group(4)

    if chain_id is None:
        chain_id = "A"
        _LOGGER.info(
            f"No chain id is provided in: {mutation_flag}. Using A as default."
        )
    if orig == target:
        _LOGGER.warning(f"equivalent mutation detected in {mutation_flag}. Making it WT.")
        return Mutation(None, "WT", None, None)

    return Mutation(orig, target, chain_id, res_idx)

def is_valid_mutation(mut: Mutation, stru: es.Structure) -> bool:
    """Checks if the supplied Mutation() namedtuple is valid according to the below criteria:
    (Non-WT cases)
    Mutation.orig: if match the original residue in the {stru}
    Mutation.target: a one-letter amino-acid code in the allowed list & different from orig
    Mutation.chain_id: a single letter, should exist in {stru}
    Mutation.res_idx: a 1-indexed int(), should exist in {stru}

    Args:
        mut: The Mutation() namedtuple to be judged.
        stru: the reference structure

    Raise:
        pass
    Returns:
        True if the Mutation() passes all checks, False if not.
    """
    #TODO
    if (not isinstance(mut.orig, str) or not isinstance(mut.target, str)
            or not isinstance(mut.chain_id, str)):
        return False

    if mut.orig not in ONE_LETTER_AA_MAPPER and mut.orig != "X":
        return False

    if mut.target not in ONE_LETTER_AA_MAPPER:
        return False

    if mut.target == mut.orig:
        return False

    if (mut.chain_id.upper() not in string.ascii_uppercase) and len(mut.chain_id.strip()):
        return False

    if not isinstance(mut.res_idx, int) or mut.res_idx < 1:
        return False

    return True




def generate_all_mutations(
    structure: es.Structure, ) -> Dict[Tuple[str, int], List[Mutation]]:
    """Creates all possible mutations for a given Structure() object. Puts all the mutations into a dict()
    where the (key, value) pairs are ((chain_id, residue), List[Mutation]). The List[Mutation] is all mutations
    from the existing residue to the other 20 residues. For a given enzyme with N residues, there will be a total
    of N*20 Mutation() objects in the dict().


    Args:
        structure: The Structure() object to build mutations from.

    Returns:
        A dict() with (key, value) pairs of ((chain_id, residue), List[Mutation]).
    """
    result = dict()
    residues: List[es.Residue] = structure.residues
    residues = list(filter(lambda rr: rr.is_canonical(), residues))
    for res in residues:
        (chain_id, num) = res.key()
        orig: str = res.name
        if len(orig) != 1:
            orig = chem.convert_to_one_letter(orig)
        result[res.key()] = list(
            map(
                lambda ch: Mutation(orig=orig, target=ch, chain_id=chain_id, res_idx=num),
                chem.one_letters_except(orig),
            ))

    # a last check
    for mut_list in result.values():
        for mt in mut_list:
            assert is_valid_mutation(mt)
    return result


def size_increase(mut: Mutation) -> bool:
    """Checks if the mutation described in the supplied namedtuple describes an increase in
        size for the selected residue. DOES NOT check if the Mutation() is valid.

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        If the described mutation leads to an increase in size.
    """

    return chem.RESIDUE_VOLUME_MAPPER[mut.target] > chem.RESIDUE_VOLUME_MAPPER[mut.orig]


def size_decrease(mut: Mutation) -> bool:
    """Checks if the mutation described in the supplied namedtuple describes an decrease in
        size for the selected residue. DOES NOT check if the Mutation() is valid.

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        If the described mutation leads to a decrease in size.
    """

    return chem.RESIDUE_VOLUME_MAPPER[mut.target] < chem.RESIDUE_VOLUME_MAPPER[mut.orig]


def polarity_change(mut: Mutation) -> bool:
    """Checks if the mutation described in the supplied namedtuple describes a change in polarity.
    DOES NOT check if the Mutation() is valid.

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        If the described mutation leads to a change in polarity.
    """
    return chem.residue.residue_polarity(mut.orig) != chem.residue.residue_polarity(
        mut.target)


def same_polarity(mut: Mutation) -> bool:
    """Checks if the mutation described in the supplied namedtuple does not describe a change in
    polarity. DOES NOT check if the Mutation() is valid.

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        If the described mutation leads to no change in polarity.
    """
    return not polarity_change(mut)

"""Submodule that defines the Mutation namedtuple() which describes a single point mutation in an 
enzyme. Additionally provides utility functions for determnining if Mutation()'s are defined
and if they satisfy certain change requirements.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-15
"""
import string
from collections import namedtuple
from typing import List, Dict, Tuple

import enzy_htp.structure as es
import enzy_htp.chemical as chem
from enzy_htp.chemical import ONE_LETTER_AA_MAPPER

Mutation = namedtuple("Mutation", "orig target chain_id res_num")
Mutation.__doc__ = f"""Named tuple representing a single point mutation in an enzyme.
   

	Attributes:
   		orig: the one-letter code of the original amino acid. Can be in [ "{", ".join(ONE_LETTER_AA_MAPPER.keys())}"] or "X".
		target: the one-letter code of the target mutation. Can be in [ "{", ".join(ONE_LETTER_AA_MAPPER.keys())}"] or "X".
		chain_id: a single capital letter.
		res_num: the 1-indexed int() of the residue to Mutate
"""


def valid_mutation(mut: Mutation) -> bool:
    """Checks if the supplied Mutation() namedtuple is valid according to the below criteria:
    Mutation.orig: a one-letter amino-acid code.
    Mutation.target: a one-letter amino-acid code different thatn Mutation.orig.
    Mutation.chain_id: a single letter, can also be blank or whitespace.
    Mutation.res_num: a 1-indexed int().

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        True if the Mutation() passes all checks, False if not.
    """
    if (
        not isinstance(mut.orig, str)
        or not isinstance(mut.target, str)
        or not isinstance(mut.chain_id, str)
    ):
        return False

    if mut.orig not in ONE_LETTER_AA_MAPPER and mut.orig != "X":
        return False

    if mut.target not in ONE_LETTER_AA_MAPPER:
        return False

    if mut.target == mut.orig:
        return False

    if (mut.chain_id.upper() not in string.ascii_uppercase) and len(
        mut.chain_id.strip()
    ):
        return False

    if not isinstance(mut.res_num, int) or mut.res_num < 1:
        return False

    return True


def generate_all_mutations(
    structure: es.Structure,
) -> Dict[Tuple[str, int], List[Mutation]]:
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
    residues: List[es.Residue] = structure.residues()
    residues = list(filter(lambda rr: rr.is_canonical(), residues))
    for res in residues:
        (chain_id, num) = res.sort_key()
        orig: str = res.name
        if len(orig) != 1:
            orig = chem.convert_to_one_letter(orig)
        result[res.sort_key()] = list(
            map(
                lambda ch: Mutation(
                    orig=orig, target=ch, chain_id=chain_id, res_num=num
                ),
                chem.one_letters_except(orig),
            )
        )

    # a last check
    for mut_list in result.values():
        for mt in mut_list:
            assert valid_mutation(mt)
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
        mut.target
    )


def same_polarity(mut: Mutation) -> bool:
    """Checks if the mutation described in the supplied namedtuple does not describe a change in
    polarity. DOES NOT check if the Mutation() is valid.

    Args:
        mut: The Mutation() namedtuple to be judged.

    Returns:
        If the described mutation leads to no change in polarity.
    """
    return not polarity_change(mut)

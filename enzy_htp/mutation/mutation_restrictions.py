"""Defines the MutationRestrictions() class which is responsible for storing rules for restricting
mutations on a given structure. This class is capable of restricting all mutations on specific
residues or entire chains and additionally restricting size increase/decreases on specific 
residues due to mutation identity as well as restricting or forcing polarity changes due to 
mutation identity. 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-19
"""
from copy import deepcopy
from typing import Dict, List, Set, Tuple, Any


from enzy_htp.core import InvalidMutationRestriction
import enzy_htp.structure as es

from .mutation import Mutation


class MutationRestrictions:
    """

    Attributes:
        mapper:
        pdb:
    """

    def __init__(self, mapper: Dict[str, Dict], pdb: str):
        """Simple constructor that initializes object with mapper and pdb. Note that mapper is deepcopied."""
        self.mapper = deepcopy(mapper)
        self.pdb = pdb

    def add_restriction(self, rkey: Tuple[str, int], ukey: str, uval: Any) -> None:
        pass

    def lock_chain(self, cname: str) -> None:
        pass

    def lock_residue(self, idx: int) -> None:
        pass

    def lock_residues(self, idxs: List[int]) -> None:
        pass

    def apply(self, muts: List[Mutation]) -> List[Mutation]:
        pass


def default_restriction_dict() -> Dict:
    """Convenience function that returns the restriction dict() for a single nucleotide. Description
    of the keys in the dict():
        + locked: Bool: if the entire residue is fully locked from mutation. Overrides all other restrictions.
        + illegal_targets: List[str]: specific residue mutation targets that are banned.
        + no_size_increase: Bool: Bans mutations that lead to an increase in residue volume. CANNOT BE TRUE WHEN 'no_size_decrease' IS TRUE.
        + no_size_decrease: Bool: Bans mutations that lead to an decrease in residue volume. CANNOT BE TRUE WHEN 'no_size_increase' IS TRUE.
        + no_polarity_change: Bool: Bans mutations that lead to changes in polarity. CANNOT BE TRUE WHEN 'force_polarity_change' IS TRUE.
        + force_polarity_change: Bool: Bans mutations that do not lead to changes in polarity. CANNOT BE TRUE WHEN 'no_polarity_change' IS TRUE.

    Returns:
        A dict() with the above attributes and all keys set to False except 'illegal_targets' which is an empty, deepcopied list.
    """
    keys: List[
        str
    ] = "locked illegal_targets no_size_increase no_size_decrease no_polarity_change force_polarity_change".split()
    result: Dict = dict(zip(keys, [False] * len(keys)))
    result["illegal_targets"] = deepcopy([])
    return deepcopy(result)


def valid_restriction_dict(rdict: Dict) -> bool:
    """Checks if the restriction dict() is valid. Checks for the below requirements:
        + has all of and only the keys: 'locked', 'illegal_targets', 'no_size_increase', 'no_size_decrease', 'no_polarity_change' and 'force_polarity_change'
        + all values are bool() except 'illegal_targets' which is a list() of str()'s
        + 'no_size_increase' and 'no_size_decrease' are not simultaneously True
        + 'no_polarity_change' and 'force_polarity_change' are not simultaneously True

    Args:
        rdict: A dict() from enzy_htp.mutation.mutation_restrictions.default_restriction_dict().

    Returns:
        Whether the supplied dict() satisfies all requirements.
    """
    target_keys: Set[str] = set(
        "locked illegal_targets no_size_increase no_size_decrease no_polarity_change force_polarity_change".split()
    )
    actual_keys: Set[str] = set(list(rdict.keys()))

    if target_keys != actual_keys:
        return False

    for (
        key
    ) in "locked no_size_increase no_size_decrease no_polarity_change force_polarity_change".split():
        if not isinstance(rdict[key], bool):
            return False

    itargets: List[str] = rdict["illegal_targets"]
    if not isinstance(itargets, list):
        return False

    for it in itargets:
        if not isinstance(it, str):
            return False

    if rdict["no_size_increase"] and rdict["no_size_decrease"]:
        return False

    if rdict["no_polarity_change"] and rdict["force_polarity_change"]:
        return False

    return True


def restriction_object(pdb: str) -> MutationRestrictions:
    """Method for creating a MutationRestrictions() object for a given .pdb structure. Recommended
    way for users to generate MutationRestrctions() objects.

    Args:
        pdb: A str() with the path to a .pdb structure.

    Returns:
        An initialized MutationRestrictions() object.
    """
    struct: es.Structure = es.structure_from_pdb(pdb)
    mapper: Dict[Tuple[str, int], Dict] = dict()
    for res in struct.residues():
        if not res.is_canonical():
            continue
        mapper[(res.chain(), res.num())] = default_restriction_dict()
    return MutationRestrictions(mapper, pdb)

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

from .mutation import (
    Mutation,
    size_increase,
    size_decrease,
    polarity_change,
    same_polarity,
)


class MutationRestrictions:
    """Class responsible for encoding restrictions to candidate mutations for a given .pdb file. Supports
    methods to hold mutation restrictions, add mutation restrictions, remove mutation restrictions and
    apply the restrictions to a set of candidate mutations, leaving behind only those mutations which
    are legal. Allowed restrictions include: + locking individual residues
        + locking all residues in a chain
        + locking specific mutation targets in residues
        + locking mutation-based volume changes in residues
        + locking mutation-based polarity changes in residues

    Attributes:
        mapper_: A dict() with (key, value) pairs of (chain_name, res_num) and dict()'s from default_restriction_dict().
        pdb_: A str() to the structure file the object was generated from.
    """

    def __init__(self, mapper: Dict[str, Dict], pdb: str):
        """Simple constructor that initializes object with mapper and pdb. Note that mapper is deepcopied."""
        self.mapper_ = deepcopy(mapper)
        self.pdb_ = pdb

    def add_restriction(self, rkey: Tuple[str, int], ukey: str,
                        uval: Any) -> None:
        """Method that adds a new restriction to a given residue specified by an rkey. Checks if the new
        restriction is valid and raises an error if not.

        Args:
            rkey: A tuple of form (chain_id, residue number) specifying the residue's restrictions to update.
            ukey: The restriction key to update. See `default_restriction_dict()` for more detail.
            uvalue: The value to update the restriction key with. See `default_restriction_dict()` for more detail.

        Raises:
            An InvalidMutationRestriction if the proposed new restriction is invalid or the proposed residue does not exist.
        """
        if rkey not in self.mapper_:
            raise InvalidMutationRestriction(
                f"{rkey} not in structure found at '{self.pdb_}'.")

        self.mapper_[rkey][ukey] = uval

        if not valid_restriction_dict(self.mapper_[rkey]):
            raise InvalidMutationRestriction(
                f"Restriction (key,value) pair ({ukey}, {uval}) is invalid")

    def lock_chain(self, cname: str) -> None:
        """Method that locks mutations in all residues within the specified chain. Note that no error is raised
        if the specified chain does not exist.

        Args:
            cname: Name of the chain to lock as a str().
        """
        for rkey in self.mapper_.keys():
            if rkey[0] == cname:
                self.lock_residue(rkey)

    def lock_residue(self, rkey: Tuple[str, int]) -> None:
        """Method that locks mutations for a specific residue in the structure. Wrapper on MutationRestrictions.add_restriction
        where rkey is the specified rkey and the (key, value) pair of (ukey, uvalue) is ('locked', True).

        Args:
            rkey: A Tuple of the form (chain_id, res_num) specifying the residue to lock.
        """
        self.add_restriction(rkey, "locked", True)

    def apply(
        self, muts: Dict[Tuple[str, int], List[Mutation]]
    ) -> Dict[Tuple[str, int], List[Mutation]]:
        """Applies the restrictions specified from the object to the supplied mutations.

        Args:
            mut: A dict() containing all possible mutations for each of the residue positions. Should be generated from enzy_htp.mutation.generate_all_mutations()

                Returns:
            A dict() with the same layout as the supplied one, but only with Mutation()'s satisfying the restrictions.
        """
        for rkey in muts.keys():

            restrict: Dict = self.mapper_[rkey]
            if restrict["locked"]:
                muts[rkey] = list()
                continue

            mlist: List[Mutation] = muts[rkey]
            if restrict["illegal_targets"]:
                mlist = list(
                    filter(
                        lambda mm: mm.target not in restrict["illegal_targets"],
                        mlist))

            if restrict["no_size_increase"]:
                mlist = list(filter(lambda mm: not size_increase(mm), mlist))

            if restrict["no_size_decrease"]:
                mlist = list(filter(lambda mm: not size_decrease(mm), mlist))

            if restrict["no_polarity_change"]:
                mlist = list(filter(lambda mm: not polarity_change(mm), mlist))

            if restrict["force_polarity_change"]:
                mlist = list(filter(lambda mm: polarity_change(mm), mlist))

            muts[rkey] = mlist

        return muts

    def pdb(self) -> str:
        """Getter for the file the MutationRestriction() was generated from."""
        return self.pdb_


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
        str] = "locked illegal_targets no_size_increase no_size_decrease no_polarity_change force_polarity_change".split(
        )
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
        "locked illegal_targets no_size_increase no_size_decrease no_polarity_change force_polarity_change"
        .split())
    actual_keys: Set[str] = set(list(rdict.keys()))

    if target_keys != actual_keys:
        return False

    for (
            key
    ) in "locked no_size_increase no_size_decrease no_polarity_change force_polarity_change".split(
    ):
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
    struct: es.Structure = es.PDBParser.get_structure(pdb)
    mapper: Dict[Tuple[str, int], Dict] = dict()
    for res in struct.residues:
        if not res.is_canonical():
            continue
        
        cname:str = ''
        if res.chain:
            cname = res.chain.name

        mapper[(cname, res.idx)] = default_restriction_dict()
    return MutationRestrictions(mapper, pdb)

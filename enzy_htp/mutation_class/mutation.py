"""Submodule that defines the Mutation namedtuple() which describes a single point mutation in an 
enzyme. Additionally provides utility functions for determnining if Mutation()'s are defined
and if they satisfy certain change requirements.
Also provide utilities for mutant, which is a list of Mutation()

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
        QZ Shao <shaoqz@icloud.com>
Date: 2022-06-15
"""
import copy
import re
from typing import List, Dict, Tuple, Union
from enzy_htp.chemical.residue import (CAA_CHARGE_MAPPER, ONE_LETTER_AA_MAPPER, THREE_LETTER_AA_MAPPER, convert_to_canonical_three_letter,
                                       convert_to_one_letter, convert_to_three_letter)

from enzy_htp.core.exception import InvalidMutationFlagSyntax, InvalidMutation
from enzy_htp.core.general import get_copy_of_deleted_dict
from enzy_htp.core.logger import _LOGGER
import enzy_htp.chemical as chem
import enzy_htp.structure as es


class Mutation:
    """representing a single point mutation in an enzyme.
	Attributes:
   		orig: the three-letter code of the original amino acid. (NCAA or canonical 3-letter name for CAA)
		target: the three-letter code of the target mutation. (NCAA or canonical 3-letter name for CAA. e.g: HIS not HIE)
		chain_id: a single capital letter.
		res_idx: the 1-indexed int() of the residue to Mutate

        *In the case of WT, the tuple is defined as (None, "WT", None, None)
    """

    SUPPORTED_MUTATION_TARGET_LIST = get_copy_of_deleted_dict(ONE_LETTER_AA_MAPPER, "U")
    """The list of EnzyHTP supported mutation target-residue list.
    add upon supporting. will do ncaa in the future
    TODO should be 3 letter list tho"""

    def __init__(self, orig, target, chain_id, res_idx):
        self.orig = orig
        self.target = target
        self.chain_id = chain_id
        self.res_idx = res_idx

    @property
    def mutation_tuple(self) -> tuple:
        """return the tuple form of the mutation information"""
        return (self.orig, self.target, self.chain_id, self.res_idx)

    # == checker ==
    def is_valid_mutation(self, stru: es.Structure) -> bool:
        """Checks if the supplied Mutation() namedtuple is valid according to the below criteria:
        (Non-WT cases)
        Mutation.orig: if match the original residue in the {stru}
        Mutation.target: a one-letter amino-acid code in the allowed list & different from orig
        Mutation.chain_id: should exist in {stru}
        Mutation.res_idx: should exist in {stru}

        Args:
            mut: The Mutation() namedtuple to be judged.
            stru: the reference structure

        Raise:
            enzy_htp.core.exception.InvalidMutation
        Returns:
            True if the Mutation() passes all checks.
        """
        # WT case
        if self == (None, "WT", None, None):
            return True

        # get data type right
        # yapf: disable
        if (not isinstance(self.orig, str) or not isinstance(self.target, str)
                or not isinstance(self.chain_id, str)
                or not isinstance(self.res_idx, int)):
            raise InvalidMutation(f"wrong data type in: {self}")
        # yapf: enable

        # Mutation.chain_id, Mutation.res_idx: should exist in {stru}, should not be empty
        if self.chain_id.strip() is "":
            raise InvalidMutation(f"empty chain_id in: {self}")
        if self.chain_id not in stru.chain_mapper:
            raise InvalidMutation(f"chain id in {self} does not exist in structure (in-stru: {stru.chain_mapper.keys()})")
        if self.res_idx not in stru[self.chain_id].residue_idxs:
            raise InvalidMutation(f"res_idx in {self} does not exist in structure (in-stru: {stru[self.chain_id].residue_idx_interval()})")

        # Mutation.orig: if match the original residue in the {stru}
        real_orig = convert_to_canonical_three_letter(stru[self.chain_id].find_residue_idx(self.res_idx).name)
        if real_orig != self.orig:
            raise InvalidMutation(f"original residue does not match in: {self} (real_orig: {real_orig})")

        # Mutation.target: a one-letter amino-acid code in the allowed list & different from orig
        if self.get_target(if_one_letter=True) not in self.SUPPORTED_MUTATION_TARGET_LIST:
            raise InvalidMutation(f"unsupported target residue in: {self}")
        if self.target == self.orig:
            raise InvalidMutation(f"equivalent mutation detected in: {self}. Should be (None, \"WT\", None, None).")

        return True

    def is_wild_type(self) -> bool:
        """check if self is wild type"""
        return self.target == "WT"

    # == property getter ==
    def get_target(self, if_one_letter: bool = False) -> str:
        """get the mutation target in 1/3-letter format"""
        if if_one_letter and (self.target in THREE_LETTER_AA_MAPPER):  # for ncaa we want to use 3 letter
            return convert_to_one_letter(self.target)
        return self.target

    def get_orig(self, if_one_letter: bool = False) -> str:
        """get the original residue in 1/3-letter format"""
        if if_one_letter and (self.orig in THREE_LETTER_AA_MAPPER):  # for ncaa we want to use 3 letter
            return convert_to_one_letter(self.orig)
        return self.orig

    def get_position_key(self) -> Tuple[str, int]:
        """get the position key in a form of (chain_id, res_idx)
        from the Mutation"""
        return (self.chain_id, self.res_idx)

    def get_charge_diff(self) -> int:
        """get the charge difference before and after the mutation."""
        return CAA_CHARGE_MAPPER[self.target] - CAA_CHARGE_MAPPER[self.orig]

    def str_wo_chain(self) -> str:
        """return a string represntation without the chain id"""
        return f"{self.get_orig(if_one_letter=True)}{self.res_idx}{self.get_target(if_one_letter=True)}"

    # == editor ==
    def changed_clone(self, **kwarg):
        """make a clone of self, values may chain"""
        result = copy.deepcopy(self)
        allowed_attr_names = ["orig", "target", "chain_id", "res_idx"]
        for k, v in kwarg.items():
            if k not in allowed_attr_names:
                raise TypeError(f"Attribute {k} does not exist.")
            setattr(result, k, v)
        return result

    # special
    def __eq__(self, target) -> bool:
        """mimik the nametuple behavior. historical reason."""
        return self.mutation_tuple == target

    def __getitem__(self, idx: int) -> Union[str, int]:
        """mimik the nametuple behavior. historical reason."""
        return self.mutation_tuple[idx]

    def __hash__(self) -> int:
        """support set()"""
        return hash(self.mutation_tuple)

    def __str__(self) -> str:
        """give a string representation of the mutation"""
        if self.target == "WT":
            return "WT"
        return f"{self.get_orig(if_one_letter=True)}{self.chain_id}{self.res_idx}{self.get_target(if_one_letter=True)}"

    def __repr__(self) -> str:
        """give a string representation of the mutation"""
        return f"({repr(self.get_orig())},{repr(self.get_target())},{repr(self.chain_id)},{repr(self.res_idx)})"


# == constructor ==
def generate_from_mutation_flag(mutation_flag: str) -> Mutation:
    """XA##Y -> ("X", "Y", "A", ##)
    WT -> (None, "WT", None, None)
    XA##X -> (None, "WT", None, None)
    *we may need to support 3-letter mutation in the future like: TYQA##TYP"""

    mutation_flag = mutation_flag.strip()
    if mutation_flag == "WT":
        return Mutation(None, "WT", None, None)
    pattern = r"([A-Z])([A-Z])?([0-9]+)([A-Z])"
    flag_match = re.match(pattern, mutation_flag)
    if flag_match is None:
        raise InvalidMutationFlagSyntax(f"{mutation_flag} doesnt match ([A-Z])([A-Z])?([0-9]+)([A-Z])")

    orig = convert_to_three_letter(flag_match.group(1))
    chain_id = flag_match.group(2)
    res_idx = int(flag_match.group(3))
    target = convert_to_three_letter(flag_match.group(4))

    if chain_id is None:
        chain_id = "A"
        _LOGGER.info(f"No chain id is provided in: {mutation_flag}. Using A as default.")
    if orig == target:
        _LOGGER.warning(f"equivalent mutation detected in {mutation_flag}. Making it WT.")
        return Mutation(None, "WT", None, None)

    return Mutation(orig, target, chain_id, res_idx)


def generate_mutation_from_target_list(position: Tuple[str, int], orig_resi: str, target_list: str) -> List[Mutation]:
    """generate a list of Mutation() objects from position and a list of target residues
    Args:
        position: the sequence position of the mutation (e.g.: A 101)
        orig_resi: the original residue on this position
        target_list: a list of target residue names
    Returns:
        a list of Mutation objects that corresponds to substitution on this position
        to each residues in the target list"""
    result = []
    for target in target_list:
        result.append(
            Mutation(orig=convert_to_canonical_three_letter(orig_resi),
                     target=convert_to_canonical_three_letter(target),
                     chain_id=position[0],
                     res_idx=position[1]))
    return result


# --Mutant-- # TODO may be a class in the future
# below utilities are for a list of Mutation()
# == checker ==
def check_repeat_mutation(mutant: List[Mutation]):
    """check if there is any repeating mutation (position-wise)
    in the mutant"""
    mutation_posi = []
    for mut in mutant:
        position_key = mut.get_position_key()
        if position_key in mutation_posi:
            return True
        mutation_posi.append(position_key)
    return False


# == editor ==
def remove_repeat_mutation(mutant: List[Mutation], keep: str = "last") -> List[Mutation]:
    """remove all mutations that have repeating position in the mutant
    and keep only the {keep} specificed one.
    Args:
        mutant: the target mutant
        keep: support only 'frist' and 'last' right now"""
    # san check
    if keep not in ["first", "last"]:
        raise ValueError("Only support 'first' or 'last' right now")
    mutation_mapper = {}
    for mut in mutant:
        position_key = mut.get_position_key()
        if position_key in mutation_mapper:
            if keep == "last":
                mutation_mapper[position_key] = mut
            if keep == "first":
                continue
        mutation_mapper[position_key] = mut
    return list(mutation_mapper.values())


# == property getter ==
def get_mutant_name_tag(mutant: List[Mutation]) -> str:
    """get a string tag for the mutant
    e.g.: '_A##B_C##D'"""
    return f"_{'_'.join(str(i) for i in mutant)}"

def get_mutant_name_str(mutant: List[Mutation]) -> str:
    """get a string representation of a mutant.
    e.g.: 'A##B C##D'"""
    return " ".join(str(i) for i in mutant)

# --Mutant Space-- # TODO may be a class in the future
# below utilities are for a list of mutants (each is a list of Mutation())
# == property getter ==
def get_involved_mutation(mutant_space: List[List[Mutation]]) -> List[Mutation]:
    """get all involved single point Mutations from a mutant space"""
    result = set()
    for mutant in mutant_space:
        for mutation in mutant:
            result.add(mutation)
    return result

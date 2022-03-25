"""
Sub module that holds various functions related to mutating residue identities. Includes 
the Mutation() namedtuple as well as utility functions.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
from typing import Union, List, Tuple
from collections import namedtuple

import enzy_htp.core as core
import enzy_htp.chemical as chem
#TODO(CJ): Function that looks for good places to mutate. 
# smaller to larger
# larger to smaller
# positive to negative
# positive to neutral 
# neutral to positive 
# neutral to negative 
# negative to positive
# negative to neutral 
Mutation = namedtuple(
    "Mutation", "orig_residue chain_index residue_index target_residue"
)
"""Namedtuple that represents a mutation that is supposed to occur in a PDB."""

def mutation_to_str(mf: Mutation) -> str:
    """Converts a mutation to a string."""
    return f"{mf.orig_residue}{mf.chain_index}{mf.residue_index}{mf.target_residue}"


def decode_mutaflags(flags : Union[str, List[str]]) -> Mutation:
    """Utility function that turns a single or list of mutflag strings into Mutation objects."""
    if isinstance(flags, str):
        flags = [flags]
    decoded = []
    
    for raw in flags:
        if len(raw) < 4:
            core._LOGGER.warning(f"The supplied flag '{raw}' is too short to be valid. Must be at least 4 characters long. Skipping mutaflag...")
            continue
        orig_res = raw[0]    
        if orig_res not in chem.ONE_LETTER_AA_MAPPER and orig_res != 'X':
            core._LOGGER.warning(f"The supplied original residue '{orig_res}' is invalid. Skipping mutaflag...")
            continue
        raw: str = raw[1:]
        chain_index: str = raw[0]
        if not chain_index.upper() or not chain_index.isalpha():
            core._LOGGER.warning(f"The supplied chain index '{chain_index}' is NOT a capitlized single letter. Maybe mal-formed.")
        raw = raw[1:]
        target_res : str = raw[-1]
        raw = raw[:-1]
        if target_res not in chem.ONE_LETTER_AA_MAPPER and target_res != 'X':
            core._LOGGER.warning(f"The supplied original residue '{target_res}' is invalid. Skipping mutaflag...")
            continue
        residue_index : str = raw
        try:
            residue_index = int(residue_index)
        except Exception as err:
            core._LOGGER.error(f"Unable to convert '{residue_index}' to an int(). Skipping mutation...")
            continue
        decoded.append( Mutation( orig_residue=orig_res, chain_index=chain_index, residue_index=residue_index, target_residue=target_res))
    
    return decoded


def get_all_combinations( curr_state : List[Tuple[str,str,int]], restrictions : List[Union[int,Tuple[int,int]]] = None, allow_U : bool = False ) -> List[Mutation]:
    """Given the current residue state derived from Structure.residue_state() and a list of index restrictions, generates all possible mutations."""
    all_combos = []
    illegal = set()
    if restrictions:
        if isinstance(restrictions,int):
            illegal.add( restrictions )
        else:
            for rr in restrictions:
                if isinstance(rr,int):
                    illegal.add(rr)
                elif len(rr) == 2:
                    illegal |= set(list(range(rr[0], rr[1]+1)))
                else:
                    core._LOGGER.error(f"Invalid restriction syntax: {restrictions}")

    for (chain, res_name, index)  in curr_state:
        if index in illegal:
            continue
        
        if len(res_name) != 1:
            core._LOGGER.warning(f"get_all_combinations() expects residue states to have one letter residue names. Continuing..")
            continue

        for other in chem.one_letters_except( res_name ):
            all_combos.append(Mutation(orig_residue=res_name,chain_index=chain,residue_index=index,target_residue=other))

    return all_combos


"""Functions related to decode the position_pattern. More specifically it should be called residue
sequence position selection pattern.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""
from typing import List
from enzy_htp.structure import Structure, Residue
from enzy_htp.structure.structure_selection import select_stru


def decode_position_pattern(stru: Structure, pattern: str, if_name: bool = False) -> List[tuple]:
    """decode pattern of residue position selection
    TODO support customized position selector
    Args:
        stru: the Structure object of reference
        pattern: a pymol-like syntax to select residue positions
                 NOTE: the result is different from pymol's result
                 that all non polypeptide part are filtered
    Returns:
        (chain_id, resi_idx) to indicate a mutation position
        ((chain_id, resi_idx), resi_name)  if_name=True"""
    selection_obj = select_stru(stru, pattern)
    result_residue: Residue = selection_obj.involved_residues
    result_residue = filter(lambda i: i.is_canonical() or i.is_modified(),
                            result_residue)  # because it is impossible and dont make sense to mutate non-polypeptide
    if if_name:
        result = [(x.key(), x.name) for x in result_residue]
    else:
        result = [x.key() for x in result_residue]

    return result

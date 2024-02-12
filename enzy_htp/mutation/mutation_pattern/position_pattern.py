"""Functions related to decode the position_pattern. More specifically it should be called residue
sequence position selection pattern.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""
import re
import numpy as np
from plum import dispatch
from typing import Iterable, List, Tuple

import enzy_htp.core.math_helper as mh
from enzy_htp.core.general import split_but_brackets
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Structure, Residue
from enzy_htp.structure.structure_selection import select_stru


def decode_position_pattern(stru: Structure, pattern: str, if_name: bool = False) -> List[tuple]:
    """decode pattern of residue position selection
    TODO support customized position selector
    Args:
        stru: the Structure object of reference
        pattern: 
            1. a pymol-like syntax to select residue positions
               NOTE: the result is different from pymol's result
               that all non polypeptide part are filtered
               example: "resi 254 around 5"
            2. a pattern starts with "$" will be recognized to use
               an internal EnzyHTP function to selected positions.
               example: "$ef_hotspot(B.254.CAE, B.254.H2, (0,10))" means residues
               within a cone region that formed by points that have angle
               <= 10 degree than vector from B.254.CAE to B.254.H2.

    Returns:
        (chain_id, resi_idx) to indicate a mutation position
        ((chain_id, resi_idx), resi_name)  if_name=True"""
    # 1. dispatch pattern type
    pattern = pattern.strip()
    if pattern.startswith("$"):
        result_residue = decode_builtin_function(stru, pattern)
    else:
        result_residue = decode_pymol_selection(stru, pattern)

    if if_name:
        result = [(x.key(), x.name) for x in result_residue]
    else:
        result = [x.key() for x in result_residue]

    return result

def decode_builtin_function(stru: Structure, pattern: str) -> Iterable[Residue]:
    """decode the position pattern in builtin function format.
    example: '$ef_hotspot('B.254.CAE', 'B.254.H2', (0,10))'"""
    pattern_format = r"\$(\w+)\((.+)+\)"
    funcname, args = re.match(pattern_format, pattern).groups()
    
    func = PATTERN_BUILTIN_FUNCTION_MAPPER[funcname]
    args = [eval(i.strip()) for i in split_but_brackets(args, ",")]

    return func(stru, *args)


def decode_pymol_selection(stru: Structure, pattern: str) -> Iterable[Residue]:
    """decode the position pattern in pymol selection format.
    example: 'resi 254 around 5'"""
    selection_obj = select_stru(stru, pattern)
    result_residue: Residue = selection_obj.involved_residues
    result_residue = filter(lambda i: i.is_canonical() or i.is_modified(),
                            result_residue)  # because it is impossible and dont make sense to mutate non-polypeptide
    return result_residue

@dispatch
def ef_hotspot(stru: Structure, atom_1: str, atom_2: str, cutoff: Tuple) -> List[Residue]:
    """find mutation sites that are most impactful to
    the internal electric field of the enzyme in terms of
    its stablization effect on a bond (defined by {atom_1, atom_2}).
    The sites are calculated to be within a cone region defined by
    <the angle between the dipole and the vector points from the
    dipole center to residue in interval {cutoff}>
    Args:
        stru:
            the target enzyme as a pre-reaction complex
        atom_1, atom_2:
            the target bond defined by 2 atoms (atom_1 -> atom_2)
        cutoff:
            the degree cutoff value for finding the mutations site.
    Returns:
        positions as a list of Residue()"""
    result = []
    # init
    low, high = cutoff
    atom_1 = np.array(stru.get(atom_1).coord)
    atom_2 = np.array(stru.get(atom_2).coord)
    p2 = (atom_1 + atom_2) / 2
    p1 = atom_2
    # collect
    for res in stru.amino_acids:
        # 1. find mass center
        p3 = np.array(res.geom_center)
        # 2. calc angle and collect
        angle = mh.get_angle(p1, p2, p3)
        if angle <= high and angle >= low:
            _LOGGER.debug(f"found {res.key_str} - {angle}°")
            result.append(res)
    return result

@dispatch
def ef_hotspot(stru: Structure, dipole_vec: Tuple, dipole_center: Tuple, cutoff: Tuple) -> List[Residue]:
    """find mutation sites that are most impactful to the internal 
    electric field of the enzyme in terms of its stablization effect
    on a dipole (defined by {dipole_vec, dipole_center}).
    The sites are calculated to be within a cone region defined by
    <the angle between the dipole and the vector points from the
    dipole center to residue in interval {cutoff}>
    Args:
        stru:
            the target enzyme as a pre-reaction complex
        dipole_vec, dipole_center:
            the target bond dipole defined by a vector and the center location.
        cutoff:
            the degree cutoff value for finding the mutations site.
    Returns:
        positions as a list of Residue()"""
    result = []
    # init
    low, high = cutoff
    p1 = np.array(dipole_center)
    d1 = np.array(dipole_vec)
    # collect
    for res in stru.amino_acids:
        # 1. find mass center
        p2 = np.array(res.geom_center)
        d2 = p2 - p1
        # 2. calc angle and collect
        angle = mh.get_angle_vec(d1, d2)
        if angle <= high and angle >= low:
            _LOGGER.debug(f"found {res.key_str} - {angle}°")
            result.append(res)
    return result
    

PATTERN_BUILTIN_FUNCTION_MAPPER = {
    "ef_hotspot" : ef_hotspot,
}


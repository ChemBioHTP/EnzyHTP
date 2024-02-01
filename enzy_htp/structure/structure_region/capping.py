"""This submodule features the following function:

    + capping_with_residue_terminals()

This module was created out of necessity to allow type hinting 
at lower levels of the enzy_htp.structure.structure_region module.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2024-1-1
"""
from __future__ import annotations
from abc import abstractmethod, ABC
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import List, Tuple, Dict

from ..structure import Structure, Residue, Atom

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.math_helper import rotation_matrix_from_vectors
from enzy_htp import chemical as chem


from .structure_region import StructureRegion
# region == APIs ==

from .residue_caps import (
    cap_residue_free_terminal
)


def capping_with_residue_terminals(raw_region:StructureRegion,
                                   nterm_cap:str=None,
                                   cterm_cap:str=None,
                                   return_copy: bool= False,
                                   **kwargs) -> Union[None, StructureRegion]:
    """cap the raw region composed by whole residues
    only. In this case, capping only on the terminal
    of the residue. (keyword: res_ter_cap)
    Args:
        raw_region: the target StructureRegion object
        nterm_cap: the cap added to the N-ter. Defaults to acetate if None.
        cterm_cap: the cap added to the C-ter. Defaults to methylamide if None.
        return_copy: whether edit in place or return an capped copy.
    Return:
        (if not return_copy, change the region in place)
        capped_region: the StructureRegion after capping."""
    # init
    if nterm_cap is None:
        nterm_cap = "COCH3"

    if cterm_cap is None:
        cterm_cap = "NHCH3"

    # san check (only allow whole residues in the region)
    if not raw_region.is_whole_residue_only():
        _LOGGER.error(f"Found incomplete residue in region. ({raw_region}) "
                      "residue terminal capping only support capping for a region"
                      " composed by whole residues only.")
        raise ValueError

    new_atoms:List[Atom] = list()
    for res in raw_region.involved_residues:
        
        if raw_region.needs_nterm_cap(res):
            new_atoms.extend(cap_residue_free_terminal(res, nterm_cap, "nterm"))

        new_atoms.extend(res.atoms)

        if raw_region.needs_cterm_cap(res):
            new_atoms.extend(cap_residue_free_terminal(res, cterm_cap, "cterm"))

    if return_copy:
        new_region = raw_region.clone()
        new_region.atoms = new_atoms 
        return new_region
    else:
        raw_region.atoms = new_atoms 

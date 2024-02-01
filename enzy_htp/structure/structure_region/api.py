"""Module API defines creation functions that generate compatible, robust StructureRegion objects from 
a variety of inputs. Current methods include:

    + create_region_from_selection_pattern()
    + create_region_from_full_stru()

This is primarily how users should create StructureRegion()'s.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-12-30
"""
from __future__ import annotations
from collections import defaultdict
from copy import deepcopy
import numpy as np
from typing import Callable, Dict, List, Tuple, Union

from ..structure import Structure, Atom
from ..residue import Residue
from ..noncanonical_base import NonCanonicalBase

from .structure_region import StructureRegion

from .residue_caps import ResidueCap

from .capping import capping_with_residue_terminals

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.math_helper import round_by, is_integer

def create_region_from_selection_pattern(
        stru: Structure,
        pattern: str,
        capping_method: str = "res_ter_cap",
        **kwargs,
    ) -> StructureRegion:
    """Create StructureRegion from selection pattern and cap terminal Atom()'s/Residue()'s.

    Args:
        stru:
            the structure that contains the target region
        pattern:
            the selection pattern in pymol style that defines the region
        capping_method:
            the keyword of the capping method.
            See CAPPING_METHOD_MAPPER for more info
        (method specific options)
        "res_ter_cap"
            nterm_cap: the name of the cap added to the N-ter
            cterm_cap: the name of the cap added to the C-ter
    Returns:
        the StructureRegion"""
    # TODO figure out the import design. We probably need to put pymol down?
    from ..structure_selection import select_stru
    # select
    stru_sele = select_stru(stru, pattern)
    raw_region = StructureRegion(atoms=stru_sele.atoms)
    # capping
    if capping_method not in CAPPING_METHOD_MAPPER:
        _LOGGER.error(f"capping method ({capping_method}) not supported. Supported: {CAPPING_METHOD_MAPPER.keys()}")
        raise ValueError
    capping_func = CAPPING_METHOD_MAPPER[capping_method]
    #TODO(CJ): add the ole logic in here about return_copy
    capping_func(raw_region, **kwargs)

    return raw_region

def create_region_from_full_stru(stru: Structure) -> StructureRegion:
    """Create StructureRegion of the full structure. Essentially acts as a constructor/translation function
    that converts Structure to StructureRegion.

    Returns:
        the StructureRegion"""
    return StructureRegion(
        atoms = stru.atoms
    )

CAPPING_METHOD_MAPPER: Dict[str, Callable[[StructureRegion], StructureRegion]] = {
    "res_ter_cap" : capping_with_residue_terminals,
    "residue_terminal_capping" : capping_with_residue_terminals,
}

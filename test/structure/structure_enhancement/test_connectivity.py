"""Testing the connectivity API in enzy_htp.structure.structure_enhancement
Author: Sebastian Stull <sebastian.l.stull@vanderbilt.edu>
Date: 2025-02-13
"""

import os
import copy
import pytest

from enzy_htp.preparation.clean import remove_solvent
import enzy_htp.structure as struct
from enzy_htp.structure.structure_enchantment import connectivity
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER
import enzy_htp.structure.structure_region as stru_regi
from enzy_htp.structure.structure_region import capping
import enzy_htp.structure.structure_selection as stru_sele


CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
NCAA_LIB = f"{CURR_DIR}/../ncaa_lib"


def test_connectivity_maa():
    test_stru = struct.PDBParser().get_structure(f"{DATA_DIR}/3FCR_modified.pdb")
    test_stru.assign_ncaa_chargespin({"LLP": (-2, 1)})
    remove_solvent(test_stru)

    connectivity.init_connectivity(test_stru)

    fs.safe_rm(f"{NCAA_LIB}/LLP_any.prepin")

    assert test_stru.modified_residue[0].is_connected()


def test_connected_structure_deepcopy():
    """Test to make sure deepcopy works on connected structures.
    NOTE: This test used to fail because the default deep copier walks the entire
    connectivity graph depth-first, so on large molecules it can exceed Python's
    recursion limit and raise RecursionError (even though cycles are memoized).
    This should be fixed now by the change in DoubleLinkedNode on 2025-08-05."""
    # Load the same structure that causes the error in integration tests
    test_stru = struct.PDBParser().get_structure(f"{DATA_DIR}/3FCR_modified.pdb")
    test_stru.assign_ncaa_chargespin({"LLP": (-2, 1), "RLP": (-2, 1)})
    remove_solvent(test_stru)
    
    # Test deepcopy BEFORE connectivity initialization
    copied_stru_before = copy.deepcopy(test_stru)
    
    # Now initialize connectivity (this creates the circular references)
    connectivity.init_connectivity(test_stru)
    
    # Test deepcopy AFTER connectivity initialization
    copied_stru_after = copy.deepcopy(test_stru)

    # TODO also make sure the copied structure is connected and connectivity is correct. Same 5 atoms for test.
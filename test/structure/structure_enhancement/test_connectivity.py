"""Testing the connectivity API in enzy_htp.structure.structure_enhancement
Author: Sebastian Stull <sebastian.l.stull@vanderbilt.edu>
Date: 2025-02-13
"""

import os
from enzy_htp.preparation.clean import remove_solvent
import enzy_htp.structure as struct
from enzy_htp.structure.structure_enchantment import connectivity
from enzy_htp.core import file_system as fs
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

    

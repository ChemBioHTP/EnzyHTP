"""Testing the connectivity API in enzy_htp.structure.structure_enhancement
Author: Sebastian Stull <sebastian.l.stull@vanderbilt.edu>
Date: 2025-02-13
"""

import os
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

    sele = stru_sele.select_stru(
    test_stru, "resi 288")
    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    cterm_cap = "OH"
    nterm_cap = "H"


    capping.capping_with_residue_terminals(test_region, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
    test_stru = test_region.convert_to_structure(cap_as_residue=False)

    connectivity.init_connectivity(test_stru)

    fs.safe_rm(f"{NCAA_LIB}/LLP_any.prepin")

    for atom in test_stru.modified_residue[0].atoms:
        print(atom.name, atom.connect)

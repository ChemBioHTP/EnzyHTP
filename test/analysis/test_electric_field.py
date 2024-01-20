"""Testing enzy_htp.analysis.electric_field.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-18"""
import os
import numpy as np

from enzy_htp.analysis import ele_field_strength_at_along
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_ele_field_strength_at():
    """as name. Use data from a real project using EnzyHTP 1.0
    as example."""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_mutant_101_254_frame_0.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    target_bond = (
        test_stru.ligands[0].find_atom_name("CAE"),
        test_stru.ligands[0].find_atom_name("H2")
    )

    result = ele_field_strength_at_along(
        test_stru, 
        *target_bond, 
        region_pattern="chain A and (not resi 101)",
    )

    assert np.isclose(result, -7.825612571555479, atol=0.0001)

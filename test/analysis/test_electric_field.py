"""Testing enzy_htp.analysis.electric_field.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-18"""
import os
import numpy as np

from enzy_htp.analysis import ele_field_strength_at_along, d_ele_field_upon_mutation_coarse
from enzy_htp import PDBParser
from enzy_htp.mutation import assign_mutant

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

def test_ele_field_strength_at_diff_atom_source():
    """as name. Use data from a real project using EnzyHTP 1.0
    as example. use atoms from the WT"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_mutant_101_254_frame_0.pdb")
    test_stru_wt = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    target_bond = (
        test_stru_wt.ligands[0].find_atom_name("CAE"),
        test_stru_wt.ligands[0].find_atom_name("H2")
    )

    result = ele_field_strength_at_along(
        test_stru, 
        *target_bond, 
        region_pattern="chain A and (not resi 101)",
    )

    assert np.isclose(result, -7.825612571555479, atol=0.0001)

def test_d_ele_field_upon_mutation_coarse():
    "as name"
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    dp_start = np.array(test_stru.find_residue_with_key(("B", 254)).find_atom_name("CAE").coord)
    dp_end = np.array(test_stru.find_residue_with_key(("B", 254)).find_atom_name("H2").coord)
    point_of_ef_measure = (dp_start + dp_end)/2
    vec_of_ef_measure = dp_end - dp_start

    result = d_ele_field_upon_mutation_coarse(
        test_stru,
        assign_mutant(test_stru, "{EA24R}")[0][0],
        point_of_ef_measure,
        vec_of_ef_measure,
        method = "ca_coord",
        unit = "MV/cm",
    )

    assert np.isclose(result, 4.410450235250402, atol=1e-6)

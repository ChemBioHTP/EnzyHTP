"""Testing the enzy_htp.structure.structure_constraint.create_constraint.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-02
"""
import pytest
import os

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
from enzy_htp.structure.structure_selection import select_stru
import enzy_htp.structure.structure_constraint as stru_cons

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

def test_create_distance_constraint():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_distance_constraint(
        "B.254.H2", "A.101.OE2", 2.4, test_stru)
    assert test_cons.params["amber"] == {
        "rs_filepath": "{mdstep_dir}/0.rs",
        "ialtd" : 0,
        "r1" : "x-0.25",
        "r2" : "x-0.05",
        "r3" : "x+0.05",
        "r4" : "x+0.25",
        "rk2": 200.0, "rk3": 200.,
    }
    assert test_cons.constraint_type == "distance_constraint"
    assert test_cons.atom_names == {"H2", "OE2"}
    assert test_cons.target_value == 2.4

def test_create_group_distance_constraint():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_group_distance_constraint(
        "resi 10 & n. C+CA+N",
        "resi 101 & n. C+CA+N",
        4.4, test_stru)

    assert len(test_cons.atoms) == 6
    assert test_cons.params["amber"] == {
        "rs_filepath": "{mdstep_dir}/0.rs",
        "ialtd" : 0,
        "r1" : "x-0.25",
        "r2" : "x-0.05",
        "r3" : "x+0.05",
        "r4" : "x+0.25",
        "rk2": 200.0, "rk3": 200.0,
    }
    assert test_cons.constraint_type == "group_distance_constraint"
    assert test_cons.atom_names == {"C", "CA", "N"}
    assert test_cons.target_value == 4.4

def test_create_angle_constraint():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_angle_constraint(
        "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru)
    assert test_cons.params["amber"] == {
        "rs_filepath": "{mdstep_dir}/0.rs",
        "ialtd" : 0,
        "r1" : "x-30.0",
        "r2" : "x-10.0",
        "r3" : "x+10.0",
        "r4" : "x+30.0",
        "rk2": 200.0, "rk3": 200.0,
    }
    assert test_cons.constraint_type == "angle_constraint"
    assert test_cons.atom_names == {"H2", "OE2", "CAE"}
    assert test_cons.target_value == 180.0

def test_create_dihedral_constraint():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_dihedral_constraint(
        "B.254.CAE", "B.254.H2", "A.101.OE2", "A.101.CA", 0.0, test_stru)
    assert test_cons.params["amber"] == {
        "rs_filepath": "{mdstep_dir}/0.rs",
        "ialtd" : 0,
        "r1" : "x-30.0",
        "r2" : "x-10.0",
        "r3" : "x+10.0",
        "r4" : "x+30.0",
        "rk2": 200.0, "rk3": 200.0,
    }
    assert test_cons.constraint_type == "dihedral_constraint"
    assert test_cons.atom_names == {"H2", "OE2", "CAE", "CA"}
    assert test_cons.target_value == 0.0

def test_merge_cartesian_freeze():
    """test a correct case"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons_list = [
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[0:10]),
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[5:12]),
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[20:22]),]
    merged_cons = stru_cons.merge_cartesian_freeze(test_cons_list)
    assert set(merged_cons.atoms) == set(test_stru.atoms[0:12]+test_stru.atoms[20:22])
    assert merged_cons.params == test_cons_list[0].params

def test_merge_cartesian_freeze_wrong(caplog):
    """test a wrong case"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons_list = [
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[0:10]),
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[5:12]),
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[20:22]),]
    test_cons_list[0].params = {}
    with EnablePropagate(_LOGGER):
        with pytest.raises(ValueError) as e:    
            merged_cons = stru_cons.merge_cartesian_freeze(test_cons_list)
            assert "Inconsistent params" in caplog.text

def test_merge_cartesian_freeze_bb_mix():
    """test a bb freeze case"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons_list = [
        stru_cons.create_backbone_freeze(test_stru),
        stru_cons.CartesianFreeze(atoms=test_stru.atoms[5:12]),]
    merged_cons = stru_cons.merge_cartesian_freeze(test_cons_list)
    assert merged_cons.constraint_type != "backbone_freeze"

def test_merge_cartesian_freeze_bb_only():
    """test a bb freeze case"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons_list = [
        stru_cons.create_backbone_freeze(stru=test_stru),]
    merged_cons = stru_cons.merge_cartesian_freeze(test_cons_list)
    assert merged_cons is test_cons_list[0]

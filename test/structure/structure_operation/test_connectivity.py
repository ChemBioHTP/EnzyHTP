"""Testing the enzy_htp.structure.structure_operation.connectivity.py
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-03-20
"""
import os

from enzy_htp import Structure
from enzy_htp import PDBParser
from enzy_htp.structure.structure_operation import init_connectivity

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/../../test_data"
sp = PDBParser()

def test_init_connectivity():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    init_connectivity(stru)

def test_init_connectivity_lv_1():
    """level 1 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    init_connectivity(test_stru)


def test_init_connectivity_lv_2():
    """level 2 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)
    ** use existing parm files for ligand"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    init_connectivity(test_stru)


def test_init_connectivity_lv_3():
    """level 3 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    init_connectivity(test_stru)


def test_init_connectivity_lv_4():
    """level 4 test of the init_connectivity.
    Test structure diversity:
    - 2 polypeptide chain
    - 2 substrate (CHN, CHONP), special net charge"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/puo_put.pdb")
    init_connectivity(test_stru)


def test_init_connectivity_lv_5():
    """level 5 test of the init_connectivity.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHONP)"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/3cfr-slp-pea_ah.pdb")
    init_connectivity(test_stru)


def test_init_connectivity_lv_6():
    """level 6 test of the init_connectivity.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHONP)"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/tyna_clean.pdb")
    init_connectivity(test_stru)



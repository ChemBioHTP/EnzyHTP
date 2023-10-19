"""Testing the enzy_htp.structure.structure_operation.connectivity.py
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-03-20
"""
import os
import logging

from enzy_htp import _LOGGER
from enzy_htp import PDBParser
from enzy_htp.structure import Atom, Structure
from enzy_htp.structure.structure_operation import init_connectivity

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/../../test_data"
STRU_DATA_DIR = f"{CURR_DIR}/../data"
sp = PDBParser()

def test_init_connectivity(): #TODO: finish these test
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

# TODO 
def test_init_connectivity_atom_missing_h(caplog):
    """get connect with missing H atoms"""

    original_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    stru = PDBParser().get_structure(f"{STRU_DATA_DIR}1NVG.pdb")
    atom1: Atom = stru["A"][0][0]  #N (Nter)
    atom2: Atom = stru["A"][1][0]  #N
    atom3: Atom = stru["A"][-1].find_atom_name("C")  #C (Cter)
    init_connectivity(atom1)
    init_connectivity(atom2)
    init_connectivity(atom3)
    assert len(atom1.connect) == 1
    assert len(atom2.connect) == 2
    assert len(atom3.connect) == 3
    assert "missing connecting atom " in caplog.text

    _LOGGER.setLevel(original_level)


def test_init_connectivity_atom_w_h():
    """get connect with no missing H atoms"""
    stru = PDBParser().get_structure(f"{STRU_DATA_DIR}1Q4T_peptide_protonated.pdb")
    atom1: Atom = stru["A"][0][0]  #N (Nter)
    atom2: Atom = stru["A"][1][0]  #N
    atom3: Atom = stru["A"][-1].find_atom_name("C")  #N (Cter)
    init_connectivity(atom1)
    init_connectivity(atom2)
    init_connectivity(atom3)
    assert len(atom1.connect) == 4
    assert len(atom2.connect) == 3
    assert len(atom3.connect) == 3



"""Testing the enzy_htp.structure.structure_operation.connectivity.py
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-03-20
"""
import os
import logging
from pathlib import Path

import enzy_htp.core.file_system as fs
from enzy_htp.structure import Atom, Structure
from enzy_htp.structure.structure_operation import init_connectivity
from enzy_htp import _LOGGER
from enzy_htp import PDBParser
from enzy_htp import config as eh_config

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/../../test_data"
STRU_DATA_DIR = f"{CURR_DIR}/../data"
sp = PDBParser()

def test_init_connectivity_lv_1():
    """level 1 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    init_connectivity(test_stru)

    # test
    sample_atom = test_stru.atoms[50]
    assert set(sample_atom.connect_atom_names) == set(["N", "HA", "CB", "C"])


def test_init_connectivity_lv_2():
    """level 2 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)
    ** use existing parm files for ligand"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"LIGAND" : (0, 1)})
    init_connectivity(test_stru, ncaa_lib=f"{CURR_DIR}/../../_interface/data/ncaa_lib/")

    # test
    sample_atom = test_stru.ligands[0].atoms[3]
    assert not os.path.exists("../ncaa_lib") # the default ncaa lib is not generated
    assert set(sample_atom.connect_atom_names) == set(["CAF", "NAL", "CAC"])


def test_init_connectivity_lv_3():
    """level 3 test of the init_connectivity.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"LIGAND" : (0, 1)})
    init_connectivity(test_stru)

    # test
    sample_atom = test_stru.ligands[0].atoms[3]
    assert set(sample_atom.connect_atom_names) == set(["CAF", "NAL", "CAC"])

    # clean
    temp_paths = [
        eh_config['system.NCAA_LIB_PATH'],
    ] + list(Path(eh_config['system.NCAA_LIB_PATH']).glob("*"))
    fs.clean_temp_file_n_dir(temp_paths)
    for p in temp_paths:
        assert not os.path.exists(p)


def test_init_connectivity_lv_4():
    """level 4 test of the init_connectivity.
    Test structure diversity:
    - 2 polypeptide chain
    - 2 substrate (CHN, CHONP), special net charge"""
    ncaa_lib = f"{CURR_DIR}/../../_interface/data/ncaa_lib_empty/"
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/puo_put.pdb")
    test_stru.assign_ncaa_chargespin({"PUT" : (1, 1)})
    init_connectivity(test_stru, ncaa_lib=ncaa_lib)

    # test
    temp_paths = list(Path(ncaa_lib).glob("PUT*prepin"))
    sample_atom = test_stru.ligands[1].atoms[3]
    assert set(sample_atom.connect_atom_names) == set(["C3", "H13", "H14", "N6"])
    # clean
    fs.clean_temp_file_n_dir(temp_paths)


def test_init_connectivity_lv_5(): # TODO finish this after finish moldesc for maa
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



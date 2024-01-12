"""Testing the functionality of the enzy_htp.structure.charge.py
Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-01
"""
import os
import pytest
import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Atom, Chain, Residue, PDBParser, Structure
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_operation.charge as chrg

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/../../test_data"
STRU_DATA_DIR = f"{CURR_DIR}/../data"
sp = PDBParser()

def test_init_charge_caa_atom():
    """as name. test _init_charge_caa_atom"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    test_atom = test_stru.atoms[100]
    chrg._init_charge_caa_atom(test_atom, "ff19sb")

    assert test_atom.has_init_charge()
    assert test_atom.charge == -0.3204

def test_init_charge_lv_1():
    """level 1 test of the init_charge.
    - single CAA"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    test_res = test_stru.residues[0]

    chrg.init_charge(test_res)

    assert test_res.has_init_charge()

def test_init_charge_region():
    """level 1 test of the init_charge.
    - CAA region"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "resi 2+3+4+5+253"
    )

    chrg.init_charge(test_stru_region)

    # test
    sample_atom = test_stru_region.atoms[10]
    assert sample_atom.charge == 0.1

def test_init_charge_lv_2():
    """level 2 test of the init_charge.
    Test structure diversity:
    - single polypeptide chain"""
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_wo_S.pdb")
    chrg.init_charge(test_stru)

    # test
    sample_atom = test_stru.atoms[50]
    assert sample_atom.charge == -0.24

def test_init_charge_lv_3():
    """level 3 test of the init_charge.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)
    ** use existing charge files for ligand"""
    raise Exception("TODO")
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"LIGAND" : (0, 1)})
    chrg.init_charge(test_stru, ncaa_lib=f"{CURR_DIR}/../../_interface/data/ncaa_lib/")

    # test
    sample_atom = test_stru.ligands[0].atoms[3]
    assert not os.path.exists("../ncaa_lib") # the default ncaa lib is not generated
    assert set(sample_atom.connect_atom_names) == set(["CAF", "NAL", "CAC"])

def test_init_charge_lv_4():
    """level 4 test of the init_charge.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)"""
    raise Exception("TODO")
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"LIGAND" : (0, 1)})
    chrg.init_charge(test_stru)

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

def test_init_charge_lv_5():
    """level 5 test of the init_charge.
    Test structure diversity:
    - 2 polypeptide chain
    - 2 substrate (CHN, CHONP), special net charge"""
    raise Exception("TODO")
    ncaa_lib = f"{CURR_DIR}/../../_interface/data/ncaa_lib_empty/"
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/puo_put.pdb")
    test_stru.assign_ncaa_chargespin({"PUT" : (1, 1)})
    chrg.init_charge(test_stru, ncaa_lib=ncaa_lib)

    # test
    temp_paths = list(Path(ncaa_lib).glob("PUT*prepin"))
    sample_atom = test_stru.ligands[1].atoms[3]
    assert set(sample_atom.connect_atom_names) == set(["C3", "H13", "H14", "N6"])
    # clean
    fs.clean_temp_file_n_dir(temp_paths)

def test_init_charge_lv_6(): # TODO finish this after finish moldesc for maa
    """level 6 test of the init_charge.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHONP)"""
    raise Exception("TODO")
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/3cfr-slp-pea_ah.pdb")
    chrg.init_charge(test_stru)

def test_init_charge_lv_7():
    """level 7 test of the init_charge.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHON)
    - 1 metal center (Cu)"""
    raise Exception("TODO")
    test_stru = PDBParser().get_structure(
        f"{DATA_DIR}/diversed_stru/tyna_clean.pdb")
    chrg.init_charge(test_stru)


"""Testing the enzy_htp._interface.pymol_interface.PyMolInterface class.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-16
"""
import os
import shutil
from pathlib import Path

from enzy_htp import interface
from enzy_htp import PDBParser
from enzy_htp import config as eh_config
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.file_system import clean_temp_file_n_dir
from enzy_htp import config as eh_config

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/data/"


def test_load_enzy_htp_stru():
    """Testing that Structure() can be loaded correctly"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    pi = interface.pymol
    test_session = pi.new_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_session, test_stru)
    assert pymol_obj_name == "tmp01"
    assert pymol_obj_name in session.cmd.get_object_list()
    assert test_stru.chain_names == session.cmd.get_chains(pymol_obj_name)
    assert len(session.cmd.get_model(pymol_obj_name).get_residues()) == test_stru.num_residues
    #temp files
    if _LOGGER.level > 10:
        assert not os.path.exists(f"{eh_config['system.SCRATCH_DIR']}/temp_pymol_interface.pdb")


def test_load_enzy_htp_stru_not_start_from_one():
    """Testing that Structure() can be loaded correctly if the structure is from
    a PDB that does not number from 1 both for atom and residue"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_session, test_stru)
    assert pymol_obj_name == "tmp01"
    assert pymol_obj_name in session.cmd.get_object_list()
    assert test_stru.chain_names == session.cmd.get_chains(pymol_obj_name)
    assert len(session.cmd.get_model(pymol_obj_name).get_residues()) == test_stru.num_residues
    pymol_resi_idxs = []
    session.cmd.iterate(f"{pymol_obj_name} & (n. CA|n. OAB) ", "pymol_resi_idxs.append(resi)", space=locals())
    assert pymol_resi_idxs == [str(x) for x in test_stru.residue_indexes]


def test_select_pymol_obj():
    """testing _select_pymol_obj works as expected"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_session, test_stru)
    test_pattern = "resi 3"
    assert test_stru["A"].find_residue_idx(3).atom_idx_list == pi.select_pymol_obj(session, test_pattern, pymol_obj_name)
    test_pattern = "resi 254 around 5"
    assert pi.select_pymol_obj(session, test_pattern, pymol_obj_name) == [
        138, 139, 140, 141, 142, 143, 144, 145, 165, 166, 168, 169, 175, 176, 177, 178, 179, 180, 181, 741, 742, 743, 744, 745, 772, 773,
        774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 1574, 1575, 1577, 1578, 1579, 1965, 1966, 1967,
        1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 2009, 2229, 2618, 2619, 2620, 2622, 2623, 2624, 2626,
        2732, 2733
    ]


def test_point_mutate():
    """testing that point mutations work as expected"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_mut_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S_mut.pdb")
    pi = interface.pymol
    pymol_obj_name, session = pi.load_enzy_htp_stru(stru=test_stru, session=pi.new_session())
    pi.point_mutate(pos_key=("A", 154), target="TRP", pymol_obj_name=pymol_obj_name, pymol_session=session)
    assert (len(pi.select_pymol_obj(pattern="resi 154", pymol_obj_name=pymol_obj_name, session=session)) == 
    (len(test_mut_stru["A"].find_residue_idx(154).atom_idx_list)))


def test_export_enzy_htp_stru():
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(stru=test_stru, session=test_session)
    test_save_stru = pi.export_enzy_htp_stru(pymol_obj_name, session)

    for new_res, old_res in zip(test_save_stru.residues, test_stru.residues):
        assert len(new_res.atoms) == len(old_res.atoms)


def test_export_pdb():
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(stru=test_stru, session=test_session)
    test_save_file = pi.export_pdb(pymol_obj_name=pymol_obj_name, pymol_session=session)
    test_save_stru = PDBParser().get_structure(test_save_file)
    clean_temp_file_n_dir([test_save_file, eh_config["system.SCRATCH_DIR"]])

    for new_res, old_res in zip(test_save_stru.residues, test_stru.residues):
        assert len(new_res.atoms) == len(old_res.atoms)

def test_fix_pymol_naming():
    """Test function works as expected"""
    atom_freq = {}

    pymol_pdb_path = f"{DATA_DIR}pymol_TRP.pdb"
    good_pdb_path = f"{DATA_DIR}Stru_TRP.pdb"

    good_stru = PDBParser().get_structure(good_pdb_path)
    pymol_stru = PDBParser().get_structure(pymol_pdb_path)

    for atom in good_stru.atoms:
        atom_freq[atom.name] = 1
    
    pi = interface.pymol
    pi.fix_pymol_naming(pymol_stru)

    for atom in pymol_stru.atoms:
        atom_freq[atom.name] -= 1

    for atom in atom_freq:
        assert atom_freq[atom] == 0

def test_get_exposed_residues():
    """Test function works as expected.
    answer generated by visual conformation in pymol"""
    pdb_path = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    stru = PDBParser().get_structure(pdb_path)
    pi = interface.pymol
    answer = [
        1,2,3,5,11,13,14,15,16,18,19,20,21,22,24,25,26,27,28,29,31,34,37,38,40,41,42,45,46,52,54,
        55,56,57,58,59,60,61,64,67,68,70,71,72,74,75,76,82,84,85,86,87,91,93,94,95,96,98,104,105,
        108,109,110,111,114,115,118,119,120,121,122,123,132,133,134,135,136,137,143,144,145,146,147,
        148,149,150,151,152,154,155,158,159,161,162,163,164,166,173,174,175,176,177,178,179,180,184,
        185,188,189,191,192,193,195,196,197,202,205,206,207,208,211,212,215,217,219,224,225,229,230,
        231,233,235,236,239,240,242,243,244,245,247,249,251,252,253
    ]


    results = pi.get_exposed_residues(stru, 0.1)
    for res_idx, ans_idx in zip(results, answer):
        assert res_idx == ans_idx

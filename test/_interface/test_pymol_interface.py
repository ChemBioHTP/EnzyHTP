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

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/data/"


def test_load_enzy_htp_stru():
    """Testing that Structure() can be loaded correctly"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    pi = interface.pymol
    test_session = pi.new_pymol_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_stru, test_session)
    assert pymol_obj_name == "enzy_htp_stru01"
    assert pymol_obj_name in session.cmd.get_object_list()
    assert test_stru.chain_names == session.cmd.get_chains(pymol_obj_name)
    assert len(
        session.cmd.get_model(pymol_obj_name).get_residues()) == test_stru.num_residues
    #temp files
    if _LOGGER.level > 10:
        assert not os.path.exists(f"{eh_config['system.SCRATCH_DIR']}/temp_pymol_interface.pdb")


def test_load_enzy_htp_stru_not_start_from_one():
    """Testing that Structure() can be loaded correctly if the structure is from
    a PDB that does not number from 1 both for atom and residue"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_pymol_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_stru, test_session)
    assert pymol_obj_name == "enzy_htp_stru01"
    assert pymol_obj_name in session.cmd.get_object_list()
    assert test_stru.chain_names == session.cmd.get_chains(pymol_obj_name)
    assert len(
        session.cmd.get_model(pymol_obj_name).get_residues()) == test_stru.num_residues
    pymol_resi_idxs = []
    session.cmd.iterate(f"{pymol_obj_name} & (n. CA|n. OAB) ",
                        "pymol_resi_idxs.append(resi)",
                        space=locals())
    assert pymol_resi_idxs == [str(x) for x in test_stru.residue_indexes]


def test_select_pymol_obj():
    """testing _select_pymol_obj works as expected"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_trun.pdb")
    pi = interface.pymol
    test_session = pi.new_pymol_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_stru, test_session)
    test_pattern = "resi 3"
    assert test_stru["A"].find_residue_idx(3).atom_idx_list == pi.select_pymol_obj(
        test_pattern, pymol_obj_name, session)
    test_pattern = "resi 254 around 5"
    assert pi.select_pymol_obj(test_pattern, pymol_obj_name, session) == [
        138, 139, 140, 141, 142, 143, 144, 145, 165, 166, 168, 169, 175, 176, 177, 178,
        179, 180, 181, 741, 742, 743, 744, 745, 772, 773, 774, 775, 776, 777, 778, 779,
        780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 1574, 1575, 1577, 1578, 1579,
        1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977,
        1978, 1979, 1980, 2009, 2229, 2618, 2619, 2620, 2622, 2623, 2624, 2626, 2732, 2733
    ]


def test_point_mutate():
    """testing that point mutations work as expected"""
    test_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_mut_stru = PDBParser().get_structure(f"{DATA_DIR}KE_07_R7_2_S_mut.pdb")
    pi = interface.pymol
    test_session = pi.new_pymol_session()
    pymol_obj_name, session = pi.load_enzy_htp_stru(test_stru, test_session)
    pi.point_mutate(("A", 154), "TRP", pymol_obj_name, session)
    assert len(pi.select_pymol_obj("resi 154", pymol_obj_name, session)) == (
           len(test_mut_stru["A"].find_residue_idx(154).atom_idx_list)
    )

"""Testing enzy_htp.analysis.dipole.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-16
"""
import glob
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_region import create_region_from_selection_pattern
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.analysis import bond_dipole
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

@pytest.mark.accre
def test_bond_dipole():
    """Test running bond_dipole calculation
    Using LMO-Multiwfn & Accre as an example set up"""
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_region = create_region_from_selection_pattern(
        test_stru, "resi 101+254")
    target_bond = (
        test_stru.ligands[0].find_atom_name("CAE"),
        test_stru.ligands[0].find_atom_name("H2")
    )
    test_ele_stru = ElectronicStructure(
        energy_0 = 0.0,
        geometry = test_region,
        mo = f"{DATA_DIR}KE_07_R7_2_S_101_254.chk",
        mo_parser = None,
        source ="gaussian16",
    )
    dipole = bond_dipole(test_ele_stru, *target_bond)

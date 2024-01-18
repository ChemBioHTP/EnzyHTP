"""Testing enzy_htp.analysis.dipole.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-16
"""
import glob
import pytest
import os
import numpy as np

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

def test_bond_dipole():
    """Test running bond_dipole calculation
    Using LMO-Multiwfn & Accre as an example set up"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_mutant_101_254_frame_0.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})

    test_region = create_region_from_selection_pattern(
        test_stru, "resi 101+254",
        nterm_cap = "H",
        cterm_cap = "H",)
    # adjust to align with the atom order in the example
    cap_h_1 = test_region.atoms.pop(-2)
    test_region.atoms.insert(1, cap_h_1)
    cap_h_2 = test_region.atoms.pop(-1)
    test_region.atoms.insert(15, cap_h_2)

    target_bond = (
        test_stru.ligands[0].find_atom_name("CAE"),
        test_stru.ligands[0].find_atom_name("H2")
    )
    test_ele_stru = ElectronicStructure(
        energy_0 = 0.0,
        geometry = test_region,
        mo = f"{DATA_DIR}KE_mutant_101_254_frame_0.fchk",
        mo_parser = None,
        source ="gaussian16",
    )
    result = bond_dipole(test_ele_stru, target_bond[0], target_bond[1],
                         work_dir=WORK_DIR)

    assert np.isclose(result[0], 0.39995, atol=0.0001)
    assert np.array_equal(result[1], np.array((0.33744, -0.21187, -0.03473)))

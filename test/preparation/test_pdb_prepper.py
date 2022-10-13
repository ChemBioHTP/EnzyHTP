"""Testing for the enzy_htp.preparation.PDBPrepper() class. This is a borderline intergration test, comparing results to previous 
enzy_htp outputs.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pytest

from enzy_htp.core import file_system as fs
from enzy_htp.preparation import PDBPrepper

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
LIGAND_DIR = f"{CURR_DIR}/work_dir/"
PDB_NAME = f"FAcD-FA-ASP.pdb"
BASE_PDB = f"{DATA_DIR}/{PDB_NAME}"
PREPPER: PDBPrepper = None

fs.safe_rmdir(WORK_DIR)


def equiv_files(fname1: str, fname2: str, width: int = None) -> bool:
    """Helper method to check if two files are exactly equivalent."""
    for l1, l2 in zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2)):
        if width:
            l1 = l1[:width]
            l2 = l2[:width]

        if l1 != l2:
            print(f"'{l1}'")
            print(f"'{l2}'")
            return False
    return True


def test_prepper_instantiation():
    """Checking that __init__() method works and setups up directory correctly. Primarily checks that the appropriate files and dirctories are create Primarily checks that the appropriate files and dirctories are createdd."""
    global PREPPER
    pdb_copy = f"{WORK_DIR}/{PDB_NAME}"
    assert not os.path.isdir(WORK_DIR)
    assert not os.path.isdir(LIGAND_DIR)
    assert not os.path.exists(pdb_copy)

    PREPPER = PDBPrepper(BASE_PDB, work_dir=WORK_DIR)

    assert os.path.isdir(WORK_DIR)
    assert os.path.exists(pdb_copy)
    assert PREPPER.current_path() == pdb_copy


def test_prepper_rm_water():
    """Making sure that the water is removed correctly."""
    global PREPPER
    no_water_pdb = "FAcD-FA-ASP_rmW.pdb"
    target_no_water = f"{DATA_DIR}/{no_water_pdb}"
    actual_no_water = f"{WORK_DIR}/{no_water_pdb}"
    assert not os.path.exists(actual_no_water)

    assert PREPPER.rm_wat() == actual_no_water
    assert PREPPER.current_path() == actual_no_water

    assert os.path.exists(target_no_water)
    assert equiv_files(target_no_water, actual_no_water)

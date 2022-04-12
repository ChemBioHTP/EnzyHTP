"""Testing for the enzy_htp.preparation.PDBPrepper() class. This is a borderline intergration test, comparing results to previous 
enzy_htp outputs.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pytest

from enzy_htp.core import InvalidPH
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


def test_prepper_get_protonation_FAcD():
    """Making sure that the get_protonation() method works for the FAcD enzyme that has a ligand.."""
    global PREPPER
    pqr_name = "FAcD-FA-ASP_rmW.pqr.pdb"
    target_pqr = f"{DATA_DIR}/FAcD-FA-ASP_rmW.pqr"
    actual_pqr = f"{WORK_DIR}/{pqr_name}"
    actual_pqr_pdb_file = f"{WORK_DIR}/FAcD-FA-ASP_rmW_aH.pdb"
    target_pqr_pdb_file = f"{DATA_DIR}/FAcD-FA-ASP_rmW_aH.pdb"

    assert not os.path.exists(actual_pqr)
    assert not os.path.exists(actual_pqr_pdb_file)
    assert PREPPER.get_protonation() == actual_pqr_pdb_file
    assert os.path.exists(actual_pqr_pdb_file)
    assert os.path.isdir(LIGAND_DIR)
    assert PREPPER.current_path() == actual_pqr_pdb_file

    assert equiv_files(target_pqr, actual_pqr)
    assert equiv_files(target_pqr_pdb_file, actual_pqr_pdb_file, 66)


def test_prepper_get_protonation_invalid_pH():
    """Making sure that the get_protonation() method throws with an invalid pH value."""

    with pytest.raises(InvalidPH) as exe:
        PREPPER.get_protonation(ph=-1.0)

    assert exe
    assert exe.type == InvalidPH

    with pytest.raises(InvalidPH) as exe:
        PREPPER.get_protonation(ph=15)

    assert exe
    assert exe.type == InvalidPH


def test_generate_mutations():
    """Making sure the generate_mutations() method works correctly."""

    PREPPER.generate_mutations(n=1)
    assert len(PREPPER.mutations) == 1
    first_run = PREPPER.mutations[0]
    PREPPER.generate_mutations(n=1)
    assert len(PREPPER.mutations) == 1
    second_run = PREPPER.mutations[0]

    assert first_run == second_run


def test_apply_mutations():
    """Checking that the PDBPrepper() can apply the mutation using tleap."""
    pass
    #assert False

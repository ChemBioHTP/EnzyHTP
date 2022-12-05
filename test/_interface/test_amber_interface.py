"""Testing the enzy_htp.molecular_mechanics.AmberInterface class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""
import os
import shutil
from pathlib import Path
from typing import Union 

from enzy_htp.core import file_system as fs
from enzy_htp import interface

MM_BASE_DIR = Path(__file__).absolute().parent
MM_DATA_DIR = f"{MM_BASE_DIR}/data/"
MINIMIZE_INPUT_1 = f"{MM_DATA_DIR}/min_1.inp"
TARGET_MINIMIZE_INPUT_1 = f"{MM_DATA_DIR}/target_min_1.inp"
MINIMIZE_INPUT_2 = f"{MM_DATA_DIR}/min_2.inp"
TARGET_MINIMIZE_INPUT_2 = f"{MM_DATA_DIR}/target_min_2.inp"

ANTECHAMBER_FNAMES="ANTECHAMBER_AC.AC ANTECHAMBER_AC.AC0 ANTECHAMBER_AM1BCC.AC ANTECHAMBER_AM1BCC_PRE.AC ANTECHAMBER_BOND_TYPE.AC ANTECHAMBER_BOND_TYPE.AC0 ANTECHAMBER_PREP.AC ANTECHAMBER_PREP.AC0 ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.in sqm.out sqm.pdb".split()

def touch( fname : Union[str, Path] ):
    if not isinstance(fname, Path):
        fpath = Path( fname )
    else:
        fpath = fname

    if not fpath.exists():
        fh = open( fpath, 'w' )
        fh.close()

def exe_exists(exe: str) -> bool:
    """Helper method that checks if an executable exists."""
    return shutil.which(exe) is not None


def files_equivalent(fname1: str, fname2: str) -> bool:
    """Helper method that checks if the contents of two files are equivalent."""
    lines1, lines2 = fs.lines_from_file(fname1), fs.lines_from_file(fname2)

    if len(lines1) != len(lines2):
        return False

    for l1, l2 in zip(lines1, lines2):
        if l1 != l2:
            return False

    return True


def test_write_minimize_input_file():
    """Testing that minimization input files are generated correctly."""
    ai = interface.amber
    assert not Path(MINIMIZE_INPUT_1).exists()
    assert not Path(MINIMIZE_INPUT_2).exists()
    ai.write_minimize_input_file(MINIMIZE_INPUT_1, 2000)
    assert files_equivalent(MINIMIZE_INPUT_1, TARGET_MINIMIZE_INPUT_1)
    ai.write_minimize_input_file(MINIMIZE_INPUT_2, 1000)
    assert files_equivalent(MINIMIZE_INPUT_2, TARGET_MINIMIZE_INPUT_2)
    fs.safe_rm(MINIMIZE_INPUT_1)
    fs.safe_rm(MINIMIZE_INPUT_2)
    assert not Path(MINIMIZE_INPUT_1).exists()
    assert not Path(MINIMIZE_INPUT_2).exists()


def test_minimize_structure():
    """Testing that a structure can be minimized with the AmberInterface.minimize_structure() method."""
    assert False 


def test_build_param_files():
    """Testing that the AmberInterface.build_param_files() method works correctly."""
    ai = interface.amber
    assert False 


#def test_build_ligand_param_files():
#    """Testing that the AmberInterface.build_ligand_param_files() method works correctly."""
#    #TODO(CJ): at present this only checks if the files are made, not if they are correct
#    lig_frcmod = f"{MM_DATA_DIR}/atp.frcmod"
#    lig_prepin = f"{MM_DATA_DIR}/atp.prepin"
#    fs.safe_rm(lig_frcmod)
#    fs.safe_rm(lig_prepin)
#    ai = interface.amber
#    ai.build_ligand_param_files([str(f"{MM_BASE_DIR}/data/atp.pdb")], [0])
#    assert os.path.exists(lig_frcmod)
#    assert os.path.exists(lig_prepin)
#    fs.safe_rm(lig_frcmod)
#    fs.safe_rm(lig_prepin)


def test_parse_fmt_uppercase():
    """Testing that the AmberInterface.parse_fmt() method works correctly for uppercase inputs."""
    ai = interface.amber
    assert ai.parse_fmt('%FORMAT(20A4)') == (str, 4)
    
    assert ai.parse_fmt('%FORMAT(5E16.8)') == (float, 16)

    assert ai.parse_fmt('%FORMAT(10I8)') == (int, 8)

def test_parse_fmt_lowercase():
    """Testing that the AmberInterface.parse_fmt() method works correctly for lowercase inputs."""
    ai = interface.amber
    assert ai.parse_fmt('%FORMAT(20a4)') == (str, 4)
    
    assert ai.parse_fmt('%FORMAT(5e16.8)') == (float, 16)

    assert ai.parse_fmt('%FORMAT(10i8)') == (int, 8)

def test_parse_fmt_bad_input():
    """Testing that the AmberInterface.parse_fmt() method returns the correct (None,-1) for bad inputs."""

    ai = interface.amber

    assert ai.parse_fmt('') == (None,-1)
    
    assert ai.parse_fmt('20a4') == (None,-1)
    
    assert ai.parse_fmt('FORMAT(20a4)') == (None,-1)
    
    assert ai.parse_fmt('%FORMAT(20a4') == (None,-1)


def test_remove_antechamber_temp_files_cwd():
    """Testing that the AmberInterface.remove_antechamber_temp_files() method works when the current working directory is used."""
    for acf in ANTECHAMBER_FNAMES:
        touch( acf )           
        assert Path(acf).exists()

    ai = interface.amber
    ai.remove_antechamber_temp_files()

    for acf in ANTECHAMBER_FNAMES:
        assert not Path(acf).exists()

def test_remove_antechamber_temp_files_non_cwd():
    """Testing that the AmberInterface.remove_antechamber_temp_files() method works when a directory other than the current working directory is used."""
    dirname = Path('ac-test-dir/' )

    assert not dirname.is_dir()
    dirname.mkdir()

    dir_ac_fnames = list(map(lambda ll: dirname / ll,  ANTECHAMBER_FNAMES))

    for acf in dir_ac_fnames:
        touch( acf )           
        assert acf.exists()

    ai = interface.amber
    ai.remove_antechamber_temp_files( str(dirname) )

    for acf in dir_ac_fnames:
        assert not acf.exists()
    
    dirname.rmdir()


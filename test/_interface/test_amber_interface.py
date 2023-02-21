"""Testing the enzy_htp.molecular_mechanics.AmberInterface class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""
import os
import shutil
import pytest
from pathlib import Path
from typing import Union 

import enzy_htp.structure as struct
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


def test_build_param_files():
    """Testing that the AmberInterface.build_param_files() method works correctly. Does not check that the 
    actual output is correct, but only that the output is made.
    """
    ai = interface.amber
    outdir = Path('demo/')
    assert not outdir.is_dir()
    test_file:str = f"{MM_BASE_DIR}/data/ff.pdb"
    ai.build_param_files(test_file, str(outdir))
    
    fnames:List[str] = ["ff_ff.pdb", "ff.inpcrd", "ff.prmtop", "leap.in", "leap.out"]
    
    for fn in fnames:
        assert (outdir / fn).exists()
    
    shutil.rmtree( outdir )

def test_build_param_files_does_not_exist():
    """Testing that the AmberInterface.build_param_files() method with throw an error when the supplied .pdb file does NOT exist."""     


    dne = Path('dne.pdb')
    work_dir = Path('work/')

    assert not dne.exists()

    ai = interface.amber
    with pytest.raises(SystemExit) as exe:
        ai.build_param_files( dne, work_dir )                    

    shutil.rmtree( work_dir ) 

    assert exe


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


def test_parse_prmtop():
    """Testing that AmberInterface.parse_prmtop() method works. It essentially loads the supplied
    prmtop file and ensures that the correct number of items are extracted for the majority of keys.
    """
    ai = interface.amber
    data = ai.parse_prmtop(f"{MM_DATA_DIR}/prmtop_1")
    
    assert data
    
    n_atoms:int = 23231
    num_bond:int = 85
    num_ang:int = 186
    n_ptra:int = 190
    n_types:int = 19
    n_bond_h:int = 21805
    n_bond_a:int = 1450
    n_thet:int = 3319 
    n_phi:int = 6138
    n_phi_h:int = 6551
    nnb:int = 43036
    n_phb:int = 1
    n_res:int = 6967
    
    assert len(data['ATOM_NAME']) == n_atoms
    assert len(data['CHARGE']) == n_atoms
    assert len(data['ATOMIC_NUMBER']) == n_atoms
    assert len(data['MASS']) == n_atoms
    assert len(data['ATOM_TYPE_INDEX']) == n_atoms
    assert len(data['NUMBER_EXCLUDED_ATOMS']) == n_atoms
    assert len(data['NONBONDED_PARM_INDEX']) == n_types**2
    assert len(data['RESIDUE_LABEL']) == n_res
    assert len(data['RESIDUE_POINTER']) == n_res
    assert len(data['BOND_FORCE_CONSTANT']) == num_bond
    assert len(data['BOND_EQUIL_VALUE']) == num_bond
    assert len(data['ANGLE_FORCE_CONSTANT']) == num_ang
    assert len(data['ANGLE_EQUIL_VALUE']) == num_ang
    assert len(data['DIHEDRAL_FORCE_CONSTANT']) == n_ptra
    assert len(data['DIHEDRAL_PERIODICITY']) == n_ptra
    assert len(data['DIHEDRAL_PHASE']) == n_ptra
    assert len(data['SCEE_SCALE_FACTOR']) == n_ptra
    assert len(data['SOLTY']) == 52
    assert len(data['LENNARD_JONES_ACOEF']) == (n_types*(n_types+1))/2
    assert len(data['LENNARD_JONES_BCOEF']) == (n_types*(n_types+1))/2
    assert len(data['BONDS_INC_HYDROGEN']) == 3*n_bond_h
    assert len(data['BONDS_WITHOUT_HYDROGEN']) == 3*n_bond_a
    assert len(data['ANGLES_INC_HYDROGEN']) == 4*n_thet
    assert len(data['DIHEDRALS_INC_HYDROGEN']) == 5*n_phi_h
    assert len(data['DIHEDRALS_WITHOUT_HYDROGEN']) == 5*n_phi
    assert len(data['EXCLUDED_ATOMS_LIST']) == nnb
    assert len(data['HBOND_ACOEF']) == n_phb
    assert len(data['HBOND_BCOEF']) == n_phb
    assert len(data['HBCUT']) == n_phb
    assert len(data['AMBER_ATOM_TYPE']) == n_atoms
    assert len(data['TREE_CHAIN_CLASSIFICATION']) == n_atoms
    assert len(data['JOIN_ARRAY']) == n_atoms
    assert len(data['IROTAT']) == n_atoms
    assert len(data['RADII']) == n_atoms


def test_parse_prmtop_no_file():
    """Checking that AmberInteface.parse_prmtop() throws an error if the listed file does not exist."""
    
    dne = Path('./dne')
   
    assert not dne.exists()

    ai = interface.amber
    with pytest.raises(SystemExit) as exe:
        ai.parse_prmtop( dne )                    

    assert exe


def test_add_charges():
    """Checking that the AmberInterface.add_charges() method works correctly."""

    ai = interface.amber
    ss = struct.PDBParser.get_structure(f"{MM_BASE_DIR}/data/1h1d.pdb")
    assert not ss.has_charges()
    ai.add_charges( ss, f"{MM_BASE_DIR}/data/1h1d_prmtop" )
    assert ss.has_charges()


def test_add_charges_bad_file():
    """Checking that the AmberInterface.add_charge() method throws an error for bad inputs."""

    ai = interface.amber
    ss = struct.PDBParser.get_structure(f"{MM_BASE_DIR}/data/1h1d.pdb")

    with pytest.raises(SystemExit) as exe:
        ai.add_charges( ss, f"{MM_BASE_DIR}/data/bad_label_prmtop" )

    assert exe


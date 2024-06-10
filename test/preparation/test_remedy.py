"""Testing the functionality of various functions found in enzy_htp.preparation.remedy.py.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-28
"""
import os
import pytest


from typing import List

import enzy_htp 
from enzy_htp.chemical import (
    SeqRes
)

from enzy_htp.preparation import (
    protonate_stru
)

from enzy_htp.preparation.remedy import (
    identify_missing_residues,
    add_missing_residues
)

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
sp = enzy_htp.PDBParser()


def test_identify_missing_residues_no_missing() -> None:
    """Checking that when supplied a PDB code with no missing residues, nothing is returned."""
    assert not len(identify_missing_residues('150D'))


def test_identify_missing_residues_errors_on_bad_code() -> None:
    """Checking that when supplied an invalid PDB, an error is thrown."""
    py_err = None
    with pytest.raises(Exception) as py_err:
        identify_missing_residues('????')

    assert py_err is not None


def test_identify_missing_residues_gets_correct() -> None:
    """Checking that when given a valid PDB code, the correct SeqRes are created."""
    actual_result:List[SeqRes] = identify_missing_residues('1BX4')
    target_result:List[SeqRes] = [
        SeqRes(1,'A',1,'MET',True,None),
        SeqRes(1,'A',2,'THR',True,None),
        SeqRes(1,'A',3,'SER',True,None),
    ]
    assert actual_result == target_result

def test_add_missing_residues_unsupported_method() -> None:
    """Checking that when given an unsupported method, the add_missing_residues() function throws an error."""
    
    infile:str=f"{DATA_DIR}/missing_residues1.pdb"
    stru = sp.get_structure( infile )
    missing_residues = identify_missing_residues('2a2c')
    
    py_err = None
    with pytest.raises(Exception) as py_err:
        add_missing_residues( stru, missing_residues, 'not_supported')

    assert py_err is not None


def test_add_missing_residues_correct_behavior_case1() -> None:
    """Checking that the add_missing_residues() function works correctly for the test case of PDB code 2a2c."""
    
    infile:str=f"{DATA_DIR}/missing_residues1.pdb"
    stru = sp.get_structure( infile )
    missing_residues = identify_missing_residues('2a2c')
    
    add_missing_residues( stru, missing_residues, 'rosetta', WORK_DIR)
    assert len(stru.chains) == 4
    assert len(stru.chains[0].residues) == 478
    assert stru.chains[0].residues[0].idx == -19 
    assert stru.chains[0].residues[-1].idx == 458 

    assert len(stru.chains[1].residues) == 1 
    assert len(stru.chains[2].residues) == 1 
    assert len(stru.chains[3].residues) == 1 

def test_add_missing_residues_correct_behavior_case2_multiple_chains() -> None:
    """Checking that the add_missing_residues() function works correctly for the test case of PDB code <> which has two amino acid chains missing Residue()'s."""
    
    infile:str=f"{DATA_DIR}/missing_residues2.pdb"
    stru = sp.get_structure( infile )
    missing_residues = identify_missing_residues('3r3v')
    
    add_missing_residues( stru, missing_residues, 'rosetta', WORK_DIR)

    assert len(stru.chains) == 3
    assert len(stru.chains[0].residues) == 306 
    assert len(stru.chains[1].residues) == 306 
    assert stru.chains[0].residues[0].idx == -1
    assert stru.chains[1].residues[0].idx == -1

    target_sequence = ''.join("""GHMPDLADLFPGFGSEWINTSSGRIFARVGGDGPPLLLLHGFPQTHVMWHRVAPKLAERFKVIVADLPGY
GWSDMPESDEQHTPYTKRAMAKQLIEAMEQLGHVHFALAGHNRGARVSYRLALDSPGRLSKLAVLDILPT
YEYWQRMNRAYALKIYHWSFLAQPAPLPENLLGGDPDFYVKAKLASWTRAGDLSAFDPRAVEHYRIAFAD
PMRRHVMCEDYRAGAYADFEHDKIDVEAGNKIPVPMLALWGASGIAQSAATPLDVWRKWASDVQGAPIES
GHFLPEEAPDQTAEALVRFFSAAPGS""".splitlines())

    seq_mapper = stru.sequence
    assert seq_mapper['A'] == target_sequence, seq_mapper['A']
    assert seq_mapper['B'] == target_sequence, seq_mapper['B']

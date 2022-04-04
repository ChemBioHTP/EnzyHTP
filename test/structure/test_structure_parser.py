"""Testing the enzy_htp.structure.structure_parser function.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-31
"""
import os
import pytest
import string

from enzy_htp.core import file_system as fs
from enzy_htp.structure import structure_parser as sp 
from enzy_htp.structure import Structure, structure_from_pdb


CURRDIR = os.path.dirname(os.path.abspath(__file__))

def test_check_valid_pdb_good_input():
    """Good input for the check_valid_pdb() helper method."""
    dummy_pdb = f'{CURRDIR}/dummy.pdb'
    assert not os.path.exists( dummy_pdb )
    fs.write_lines( dummy_pdb, ['line1', 'line2'])
    assert not sp.check_valid_pdb( dummy_pdb )
    fs.safe_rm(dummy_pdb)
    assert not os.path.exists( dummy_pdb )


def test_check_valid_pdb_bad_input():
    """Good input for the check_valid_pdb() helper method."""
    # non pdb file
    txt_file = 'not_pdb.txt'
    assert not os.path.exists( txt_file )
    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb( txt_file )
    
    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1
    
	# pdb that doesn't exist
    pdb_imaginary = 'not_real.pdb'
    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb( pdb_imaginary )
    
    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

    non_ascii_pdb = f"{CURRDIR}/bad_pdb.pdb"
    assert not os.path.exists( non_ascii_pdb )
    fs.write_lines( non_ascii_pdb, ['日本人 中國的'] )
    assert os.path.exists( non_ascii_pdb )
    
    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb( non_ascii_pdb )
    
    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1
    
    fs.safe_rm( non_ascii_pdb )
    assert not os.path.exists( non_ascii_pdb )


def test_legal_chain_names():
    """Making sure all legal chains are created given a starting Residue() mapper."""
    ALL_NAMES = list(string.ascii_uppercase) 
    result1 = sp.legal_chain_names(dict())
    assert set(result1) == set(ALL_NAMES)
    
    dummy_mapper = dict(zip(ALL_NAMES,range(len(ALL_NAMES)))) 
    result2 = sp.legal_chain_names( dummy_mapper )
    assert not result2

    ALL_NAMES.remove('A')
    dummy_mapper = dict(zip(ALL_NAMES,range(len(ALL_NAMES)))) 
    result3 = sp.legal_chain_names( dummy_mapper )
    assert result3 == ['A'] 


def test_name_chains():
    """Ensuring that the name_chains() correctly names new chains."""
    assert False

def test_categorize_residues():
    assert False

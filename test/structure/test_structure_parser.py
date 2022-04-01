"""Testing the enzy_htp.structure.structure_parser function.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-31
"""
import os
import pytest

from enzy_htp.core import file_system as fs
from enzy_htp.structure import structure_parser as sp 
from enzy_htp.structure import Structure, structure_from_pdb


CURRDIR = os.path.dirname(os.path.abspath(__file__))

def test___check_valid_pdb_good_input():
    """Good input for the __check_valid_pdb() helper method."""
    dummy_pdb = f'{CURRDIR}/dummy.pdb'
    assert not os.path.exists( dummy_pdb )
    fs.write_lines( dummy_pdb, ['line1', 'line2'])
    assert not sp.__check_valid_pdb( dummy_pdb )
    fs.safe_rm(dummy_pdb)
    assert not os.path.exists( dummy_pdb )


def test___check_valid_pdb_bad_input():
    """Good input for the __check_valid_pdb() helper method."""
    # non pdb file
    txt_file = 'not_pdb.txt'
    assert not os.path.exists( txt_file )
    with pytest.raises(SystemExit) as exe:
        sp.__check_valid_pdb( txt_file )
    
    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1
    
	# pdb that doesn't exist
    pdb_imaginary = 'not_real.pdb'
    with pytest.raises(SystemExit) as exe:
        sp.__check_valid_pdb( pdb_imaginary )
    
    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1



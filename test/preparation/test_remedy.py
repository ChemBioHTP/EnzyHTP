"""Testing the functionality of various functions found in enzy_htp.preparation.remedy.py.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-28
"""
import pytest


from typing import List

import enzy_htp 
from enzy_htp.chemical import (
    SeqRes
)
from enzy_htp.preparation.remedy import (
    identify_missing_residues
)



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


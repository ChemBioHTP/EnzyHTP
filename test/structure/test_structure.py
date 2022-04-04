"""Testing the Structure() class. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
import os
import pytest
from copy import deepcopy
from enzy_htp.structure import (Structure, structure_from_pdb, Residue, Chain)



CURR_DIR = os.path.dirname(os.path.abspath( __file__ ))
TEST_DIR = os.path.dirname( CURR_DIR )


def test_load_structure():
    """"""
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    structure_from_pdb(TEST_FILE)
    #assert False

def test_structure_ctor_bad_input():
    """Testing that the Structure() ctor fails when given duplicate chains."""
    chain1 = Chain('A',list())
    chain2 = Chain('B',list())
    chain3 = Chain('A',list())

    with pytest.raises(SystemExit) as exe:
        struct = Structure([chain1,chain2,chain3])

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

def test_residue_state():
    """Checking that the Structure.residue_state() method properly returns the residue state for the Structure()."""
    start = [("A","ARG",1),("A","HIS",2),("A","ARG",3),("A","HIS",4)]
    answer = [("A","R",1),("A","H",2),("A","R",3),("A","H",4)]
    keys = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start))
    chain1 = Chain('A',list(map(lambda kk: Residue(kk, list()), keys)))
    ss = Structure([chain1]) #TODO(CJ): need error in Structure() ctor if NOT a list of chains
    print(ss.residue_state(),'here')
    assert ss.residue_state() == answer 

def test_structure_same_sequence():
    """Ensuring that the Chain.same_sequence() method works."""
    start1 = [("A","ARG",1),("A","HIS",2),("A","ARG",3),("A","HIS",4)]
    start2 = [("B","ARG",1),("B","HIS",2),("B","ARG",3),("B","HIS",4)]
    start3 = [("A","ARG",1),("A","ARG",2),("A","ARG",3),("A","HIS",4)]
    start4 = [("A","ARG",1),("A","HIS",2),("A","ARG",4),("A","HIS",4)]
    keys1 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start1))
    keys2 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start2))
    keys3 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start3))
    keys4 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start4))
    # testing the same chain
    chain1 = Chain('A',list(map(lambda kk: Residue(kk, list()), keys1)))
    chain2 = Chain('A',list(map(lambda kk: Residue(kk, list()), keys1)))
    assert chain1.same_sequence(chain2)
    assert chain2.same_sequence(chain1)
    # testing different chain name
    chain2 = Chain('B',list(map(lambda kk: Residue(kk, list()), keys2)))
    assert not chain1.same_sequence(chain2)
    assert not chain2.same_sequence(chain1)
    # testing different residue name 
    chain2 = Chain('A',list(map(lambda kk: Residue(kk, list()), keys3)))
    assert not chain1.same_sequence(chain2)
    assert not chain2.same_sequence(chain1)
    # testing different residue number 
    chain2 = Chain('A',list(map(lambda kk: Residue(kk, list()), keys4)))
    assert not chain1.same_sequence(chain2)
    assert not chain2.same_sequence(chain1)



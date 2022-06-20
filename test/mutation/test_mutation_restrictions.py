"""Testing the enzy_htp.mutation.mutation_restrictions.py submodule. This includes the MutationRestrictions() class
as well as the associated utility functions.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-20
"""
import os
from typing import List, Dict

import enzy_htp.mutation as mut
import enzy_htp.structure as es


DATA_DIR=f"{os.path.dirname(os.path.abspath(__file__))}/data/"


def test_default_restriction_dict():
    """Testing that the default dict() returned from enzy_htp.mutation.mutation_restrictions return what it should."""
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    assert isinstance(test_dict,dict)
    assert test_dict['locked'] == False
    assert isinstance(test_dict['illegal_targets'],list) 
    assert test_dict['no_size_increase'] == False
    assert test_dict['no_size_decrease'] == False
    assert test_dict['no_polarity_change'] == False
    assert test_dict['force_polarity_change'] == False

def test_default_restriction_dict_deepcopy():
    """Testing that the default dict()'s returned from enzy_htp.mutation.mutation_restrictions are deep-copied."""
    test_dict1:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict2:Dict = mut.mutation_restrictions.default_restriction_dict()
    assert id(test_dict1) != id(test_dict2)
    assert id(test_dict1['illegal_targets']) != id(test_dict2['illegal_targets'])


def test_valid_restriction_dict_returns_true():
    """Testing that the valid_restriction_dict() method returns False when it should."""
    assert mut.valid_restriction_dict(mut.mutation_restrictions.default_restriction_dict())
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = "A B C D E".split()
    assert mut.valid_restriction_dict(test_dict)

def test_valid_restriction_dict_returns_false():
    """Testing that the valid_restriction_dict() returns False when it should."""
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    del test_dict['illegal_targets']
    assert not mut.valid_restriction_dict(test_dict)
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['Aa'] = 10
    assert not mut.valid_restriction_dict(test_dict)
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['no_size_increase'] = True
    test_dict['no_size_decrease'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['no_polarity_change'] = True
    test_dict['force_polarity_change'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict:Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = [True]
    assert not mut.valid_restriction_dict(test_dict)



def test_restriction_object():
    """Testing that the restriction_object() method can create a simple MutationRestriction() object."""
    ONE_RES:str=f"{DATA_DIR}/one_res.pdb"
    mr:mut.MutationRestrictions=mut.restriction_object(ONE_RES)
    assert isinstance(mr, mut.MutationRestrictions)
    assert len(mr.mapper) == 1 
    assert mr.pdb ==  ONE_RES

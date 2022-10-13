"""Testing the enzy_htp.mutation.mutation_restrictions.py submodule. This includes the MutationRestrictions() class
as well as the associated utility functions.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-20
"""
import os
import pytest
from typing import List, Dict, Tuple

import enzy_htp.core as core
import enzy_htp.mutation as mut
import enzy_htp.structure as es

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"


def dummy_mr_object(rkeys: List[Tuple[str, int]]) -> mut.MutationRestrictions:
    """Helper method that uses a list() of rkeys tuple with (chain_id, r_index) to generate a dummy MutationRestrictions() class."""
    mapper = dict()
    for rr in rkeys:
        mapper[rr] = mut.mutation_restrictions.default_restriction_dict()

    return mut.MutationRestrictions(mapper, "dummy.pdb")


def test_default_restriction_dict():
    """Testing that the default dict() returned from enzy_htp.mutation.mutation_restrictions return what it should."""
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    assert isinstance(test_dict, dict)
    assert test_dict['locked'] == False
    assert isinstance(test_dict['illegal_targets'], list)
    assert test_dict['no_size_increase'] == False
    assert test_dict['no_size_decrease'] == False
    assert test_dict['no_polarity_change'] == False
    assert test_dict['force_polarity_change'] == False


def test_default_restriction_dict_deepcopy():
    """Testing that the default dict()'s returned from enzy_htp.mutation.mutation_restrictions are deep-copied."""
    test_dict1: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict2: Dict = mut.mutation_restrictions.default_restriction_dict()
    assert id(test_dict1) != id(test_dict2)
    assert id(test_dict1['illegal_targets']) != id(
        test_dict2['illegal_targets'])


def test_valid_restriction_dict_returns_true():
    """Testing that the valid_restriction_dict() method returns False when it should."""
    assert mut.valid_restriction_dict(
        mut.mutation_restrictions.default_restriction_dict())
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = "A B C D E".split()
    assert mut.valid_restriction_dict(test_dict)


def test_valid_restriction_dict_returns_false():
    """Testing that the valid_restriction_dict() returns False when it should."""
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    del test_dict['illegal_targets']
    assert not mut.valid_restriction_dict(test_dict)
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['Aa'] = 10
    assert not mut.valid_restriction_dict(test_dict)
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['no_size_increase'] = True
    test_dict['no_size_decrease'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['no_polarity_change'] = True
    test_dict['force_polarity_change'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = True
    assert not mut.valid_restriction_dict(test_dict)
    test_dict: Dict = mut.mutation_restrictions.default_restriction_dict()
    test_dict['illegal_targets'] = [True]
    assert not mut.valid_restriction_dict(test_dict)


def test_restriction_object():
    """Testing that the restriction_object() method can create a simple MutationRestriction() object."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    mr: mut.MutationRestrictions = mut.restriction_object(ONE_RES)
    assert isinstance(mr, mut.MutationRestrictions)
    assert len(mr.mapper_) == 1
    assert mr.pdb() == ONE_RES


def test_mutation_restrictions_add_restriction():
    """Testing that the MutationRestrictions.add_restriction() method works for valid input cases."""
    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])
    assert not mr.mapper_[('A', 1)]['locked']
    mr.add_restriction(('A', 1), 'locked', True)
    assert mr.mapper_[('A', 1)]['locked']
    assert not mr.mapper_[('A', 1)]['illegal_targets']
    mr.add_restriction(('A', 1), 'illegal_targets', ['ASP', 'TYR'])
    assert mr.mapper_[('A', 1)]['illegal_targets'] == ['ASP', 'TYR']


def test_mutation_restrictions_add_restriction_throws():
    """Testing that the MutationRestrictions.add_restriction() throws for missing keys and conflicting flag settings.."""

    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])
    with pytest.raises(core.InvalidMutationRestriction) as exe:
        mr.add_restriction(('A', 4), 'locked', True)

    assert exe
    assert exe.type == core.InvalidMutationRestriction

    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])
    mr.add_restriction(('A', 1), 'no_size_increase', True)
    with pytest.raises(core.InvalidMutationRestriction) as exe:
        mr.add_restriction(('A', 1), 'no_size_decrease', True)

    assert exe
    assert exe.type == core.InvalidMutationRestriction


def test_mutation_restrictions_lock_chain():
    """Testing that the MuationRestrictions.lock_chain() method works correctly."""
    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])
    mr.lock_chain('A')
    assert mr.mapper_[('A', 1)]['locked']
    assert mr.mapper_[('A', 2)]['locked']
    assert mr.mapper_[('A', 3)]['locked']

    assert not mr.mapper_[('B', 1)]['locked']
    assert not mr.mapper_[('B', 2)]['locked']
    assert not mr.mapper_[('B', 3)]['locked']

    mr.lock_chain('B')
    assert mr.mapper_[('B', 1)]['locked']
    assert mr.mapper_[('B', 2)]['locked']
    assert mr.mapper_[('B', 3)]['locked']


def test_mutation_restrictions_lock_residue():
    """Testing that the MuationRestrictions.lock_residue() method works correctly."""

    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])

    assert not mr.mapper_[('A', 1)]['locked']
    mr.lock_residue(('A', 1))
    assert mr.mapper_[('A', 1)]['locked']


def test_apply():
    """Testing that the MutationRestrictions.lock_residue() method works correctly."""
    mr: mut.MutationRestrictions = dummy_mr_object([('A', 1), ('A', 2),
                                                    ('A', 3), ('B', 1),
                                                    ('B', 2), ('B', 3)])
    mutations = dict()
    mutations[('A', 1)] = list(range(10))
    mutations[('A', 2)] = list(range(10))
    mutations[('A', 3)] = list(range(10))
    mutations[('B', 1)] = list(range(10))
    mutations[('B', 2)] = [
        mut.Mutation(orig='A', target='C', chain_id='B', res_num=10),
        mut.Mutation(orig='C', target='A', chain_id='B', res_num=10)
    ]
    mutations[('B', 3)] = [
        mut.Mutation(orig='R', target='H', chain_id='B', res_num=10),
        mut.Mutation(orig='R', target='A', chain_id='B', res_num=10),
        mut.Mutation(orig='A', target='R', chain_id='B', res_num=10)
    ]
    mr.lock_chain('A')
    mr.lock_residue(('B', 1))
    mr.add_restriction(('B', 2), 'no_size_increase', True)
    mr.add_restriction(('B', 3), 'no_polarity_change', True)
    mutations = mr.apply(mutations)

    assert not mutations[('A', 1)]
    assert not mutations[('A', 2)]
    assert not mutations[('A', 3)]
    assert not mutations[('B', 1)]
    assert mutations[('B', 2)] == [
        mut.Mutation(orig='C', target='A', chain_id='B', res_num=10)
    ]
    assert mutations[('B', 3)] == [
        mut.Mutation(orig='R', target='H', chain_id='B', res_num=10)
    ]

"""Testing the Mutation() namedtuple and associated compatability functions in the enzy_htp.mutation.mutation
submodule.



Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-15
"""
from typing import List

import enzy_htp.mutation as mut



def test_valid_mutation_passes():
    """Testing cases that should pass for enzy_htp.mutation.valid_mutation."""
    mutations : List[mut.Mutation] = [
		   mut.Mutation(orig='R', target='D', chain_id='A', res_num=1),
		   mut.Mutation(orig='X', target='D', chain_id='A', res_num=1),
		   mut.Mutation(orig='R', target='D', chain_id='', res_num=1000),
		   mut.Mutation(orig='R', target='D', chain_id=' ', res_num=1000),
		   mut.Mutation(orig='U', target='W', chain_id='A', res_num=1),
        ]

    for mm in mutations:
        assert mut.valid_mutation(mm)		


def test_valid_mutation_fails():
    """Testing cases that should fail for enzy_htp.mutation.valid_mutation."""
    mutations : List[mut.Mutation] = [
		   mut.Mutation(orig='R', target='R', chain_id='A', res_num=1),
		   mut.Mutation(orig='R', target='D', chain_id='A', res_num=0),
		   mut.Mutation(orig='R', target='D', chain_id='A', res_num=0.0),
		   mut.Mutation(orig='R', target='D', chain_id=3, res_num=1000),
		   mut.Mutation(orig='R', target='D', chain_id='3', res_num=1000),
		   mut.Mutation(orig='R', target='Z', chain_id='A', res_num=1),
		   mut.Mutation(orig='Z', target='R', chain_id='A', res_num=1),
		   mut.Mutation(orig=1, target='D', chain_id='A', res_num=0),
		   mut.Mutation(orig='R', target=1, chain_id='A', res_num=0),
		   mut.Mutation(orig='R', target='D', chain_id=1, res_num=0),
		   mut.Mutation(orig='R', target='D', chain_id='A', res_num=False),
        ]

    for mm in mutations:
        assert not mut.valid_mutation(mm)		


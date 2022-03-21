"""Testing the functionality present in enzy_htp.preparation.mutate.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""

from enzy_htp.preparation import mutate as mm 

def test_mutation_attributes():
    """Checks that the Mutation namedtuple has the proper attributes"""
    mut = mm.Mutation(orig_residue='A', chain_index='B', residue_index=2, target_residue='Z')
    assert mut.orig_residue == 'A'
    assert mut.chain_index == 'B'
    assert mut.residue_index == 2
    assert mut.target_residue == 'Z'

def test_mutation_to_str():
    """Ensuring that the mutation_to_str() method works correctly."""
    mut = mm.Mutation(orig_residue='A', chain_index='B', residue_index=2, target_residue='Z')
    mut_str = mm.mutation_to_str( mut )
    assert mut_str == "AB2Z"


def test_decode_mutaflags_good_input():
    """Testing for the decode_mutaflags() method for well formed inputs."""
    result = mm.decode_mutaflags( 'XA12H'  )
    mut = mm.Mutation(orig_residue='X', chain_index='A', residue_index=12, target_residue='H')
    assert result == [mut]
    result = mm.decode_mutaflags( ['XA12H']  )
    assert result == [mut]

def test_decode_mutaflags_bad_input():
    """Testing for the decode_mutaflags() method for poorly formed inputs."""
    assert not mm.decode_mutaflags( 'XAFH'  )
    assert not mm.decode_mutaflags( ['XAFH']  )
    assert not mm.decode_mutaflags( ['XAF3']  )
    assert not mm.decode_mutaflags( ['3AFH']  )
    assert not mm.decode_mutaflags( ['31F3']  )
    assert not mm.decode_mutaflags( ['']  )

def test_get_all_combinations():
    """Testing the ability of get_all_combinations() to generate mutation combinations."""
    res_state = [('A','H', 1)]
    combos = mm.get_all_combinations(res_state)
    assert len(combos) == 20

def test_get_all_combinations_restrictions():
    """Testing the ability of get_all_combinations() to generate mutation combinations with restrictions."""
    res_state = [('A','H', 1)]

    combos = mm.get_all_combinations(res_state, [1])
    assert not combos
    
    large_res_state = list(map(lambda idx: ('A', 'H', idx), range(1,21)))
    combos = mm.get_all_combinations(large_res_state, [1,2]) 
    indices = set(list(map( lambda mut: mut.residue_index, combos)))
    assert 1 not in indices
    assert 2 not in indices

    combos = mm.get_all_combinations(large_res_state, [(1,5)]) 
    indices = set(list(map( lambda mut: mut.residue_index, combos)))
    assert 1 not in indices
    assert 2 not in indices
    assert 3 not in indices
    assert 4 not in indices
    assert 5 not in indices


    combos = mm.get_all_combinations(large_res_state, [(1,5), (10,15), 18]) 
    indices = set(list(map( lambda mut: mut.residue_index, combos)))
    assert 1 not in indices
    assert 2 not in indices
    assert 3 not in indices
    assert 4 not in indices
    assert 5 not in indices
    assert 10 not in indices
    assert 11 not in indices
    assert 12 not in indices
    assert 13 not in indices
    assert 14 not in indices
    assert 15 not in indices
    assert 18 not in indices



"""Testing the enzy_htp.mutation.general.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
        Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-01-23
"""

import os
from typing import List

import pytest

from enzy_htp.core import file_system as fs
from enzy_htp.core.exception import UnsupportedMethod
import enzy_htp.mutation.general as mg
from enzy_htp.mutation.mutation import Mutation

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"


def test_sync_mutation_over_chains():
    """test function works as expected"""
    test_mutants = [
        [Mutation(orig='ARG',target='ALA',chain_id='A',res_idx=3),
         Mutation(orig='ARG',target='TRP',chain_id='A',res_idx=4)],
        [Mutation(orig='TRP',target='GLY',chain_id='C',res_idx=1),
         Mutation(orig='TRP',target='HIS',chain_id='C',res_idx=2)]]
    test_chain_sync_list = [("A","B"), ("C","D")]
    test_chain_index_mapper = {
        "A": 0,
        "B": 10,
        "C": 20,
        "D": 100
    }
    result = mg.sync_mutation_over_chains(
        test_mutants,
        test_chain_sync_list,
        test_chain_index_mapper)
    assert len(result) == 2
    assert len(result[0]) == 4
    assert Mutation(orig='ARG', target='TRP', chain_id='B', res_idx=14) in result[0]
    assert len(result[1]) == 4
    assert Mutation(orig='TRP', target='HIS', chain_id='D', res_idx=82) in result[1]

# == TODO ==
def test_get_mutations_random_state_works():
    """Checking that the random_state function parameter is effective for controlling output in the get_mutations() method."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    assert mut.engine.get_mutations(ONE_RES, 1, list(), 100,
                                    None) == mut.engine.get_mutations(
                                        ONE_RES, 1, list(), 100, None)
    assert mut.engine.get_mutations(ONE_RES, 1, list(), 99,
                                    None) != mut.engine.get_mutations(
                                        ONE_RES, 1, list(), 100, None)


def test_get_mutations_raises():
    """Checking that the get_mutations() method throws when there are not enough valid mutations available."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"

    with pytest.raises(Exception) as exe:
        mut.engine.get_mutations(ONE_RES, 2, list(), 100, None)
    assert exe

    restrict: mut.MutationRestrictions = mut.restriction_object(ONE_RES)
    restrict.lock_residue(('A', 1))
    with pytest.raises(Exception) as exe:
        mut.engine.get_mutations(ONE_RES, 1, list(), 100, restrict)
    assert exe


def test_mutated_name():
    """Method that tests that mutated_name() is working correctly."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    target = "dne/one_res_P1G.pdb"
    actual = mut.engine.mutated_name(
        ONE_RES, "dne", [mut.Mutation(orig='P', target='G', chain_id='A', res_idx=1)])
    assert actual == target


def test_mutated_name_can_deduce():
    """Method that tests that mutated_name() is working correctly and can specifically deduce when no original residue code is given."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    target = "dne/one_res_P1G.pdb"
    actual = mut.engine.mutated_name(
        ONE_RES, "dne", [mut.Mutation(orig='X', target='G', chain_id='A', res_idx=1)])
    assert actual == target


def test_mutate_pdb_one_letter_tleap():
    """Checking that mutate_pdb() works on a pdb structure with a single letter specifically for tleap's mutation.
    Note to developers: This should be replicated for each and every new package added."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    target = f"{DATA_DIR}/one_res_P1K.pdb"

    fs.safe_rm(target)

    assert not os.path.exists(target)
    mutated = mut.mutate_pdb(ONE_RES, 1, list(), None, "tleap", None, 100)
    assert os.path.exists(target)

    lines: List[str] = fs.lines_from_file(target)
    lines = list(filter(lambda ll: ll.startswith('ATOM'), lines))

    for ll in lines:
        assert 'LYS' in ll

    fs.safe_rm(target)
    assert not os.path.exists(target)


def test_mutate_pdb_specified_mutation_tleap():
    """Checking that mutate_pdb() works on a pdb structure with a single letter specifically for tleap's mutation 
    when that mutation has been manually specified. Note to developers: This should be replicated for each and 
    every new package added."""

    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    target = f"{DATA_DIR}/one_res_P1G.pdb"
    fs.safe_rm(target)
    assert not os.path.exists(target)
    mutated = mut.mutate_pdb(
        ONE_RES, 1, [mut.Mutation(orig='P', target='G', chain_id='A', res_idx=1)], None,
        "tleap", None, 100)
    assert os.path.exists(target)

    lines: List[str] = fs.lines_from_file(target)
    lines = list(filter(lambda ll: ll.startswith('ATOM'), lines))

    for ll in lines:
        assert 'GLY' in ll

    fs.safe_rm(target)
    assert not os.path.exists(target)


def test_mutate_pdb_raises_unsupported_method_exception():
    """Testing that the mutate_pdb() method will raise an UnsupportedMethod() Excpetion when given a nonsense 
    value for the engine."""

    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    with pytest.raises(UnsupportedMethod) as exe:
        mut.mutate_pdb(ONE_RES, 1, list(), None, "doesnt-exist")

    assert exe
    assert exe.type == UnsupportedMethod

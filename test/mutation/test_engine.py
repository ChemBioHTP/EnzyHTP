"""Testing the enzy_htp.mutation.mutate.poy submodule.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-16
"""

import os
from typing import List

import pytest

import enzy_htp.mutation as mut
import enzy_htp.structure as es
from enzy_htp.core import file_system as fs
from enzy_htp.core.exception import UnsupportedMethod

from enzy_htp import interface

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"


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
        ONE_RES, "dne", [mut.Mutation(orig='P', target='G', chain_id='A', res_num=1)])
    assert actual == target


def test_mutated_name_can_deduce():
    """Method that tests that mutated_name() is working correctly and can specifically deduce when no original residue code is given."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    target = "dne/one_res_P1G.pdb"
    actual = mut.engine.mutated_name(
        ONE_RES, "dne", [mut.Mutation(orig='X', target='G', chain_id='A', res_num=1)])
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
        ONE_RES, 1, [mut.Mutation(orig='P', target='G', chain_id='A', res_num=1)], None,
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

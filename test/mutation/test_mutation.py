"""Testing the Mutation() namedtuple and associated compatability functions in the enzy_htp.mutation.mutation
submodule.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
        QZ Shao <shaoqz@icloud.com>
Date: 2022-06-15
"""
import logging
import os
import pytest
from typing import List

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InvalidMutation
import enzy_htp.structure as es
import enzy_htp.mutation as mut

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"


def test_generate_from_mutation_flag_wt(caplog):
    """test if the function behave as expected on WT"""
    assert mut.generate_from_mutation_flag("WT") == mut.Mutation(orig=None,
                                                                 target='WT',
                                                                 chain_id=None,
                                                                 res_idx=None)

    original_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    assert mut.generate_from_mutation_flag("RA154R") == mut.Mutation(orig=None,
                                                                     target='WT',
                                                                     chain_id=None,
                                                                     res_idx=None)
    assert "equivalent mutation detected" in caplog.text

    _LOGGER.setLevel(original_level)


def test_generate_from_mutation_flag_default_chainid(caplog):
    """test if the function behave as expected on default chain id"""
    original_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    assert mut.generate_from_mutation_flag("R154A") == mut.Mutation(orig='ARG',
                                                                    target='ALA',
                                                                    chain_id='A',
                                                                    res_idx=154)
    assert " Using A as default." in caplog.text

    _LOGGER.setLevel(original_level)


def test_is_valid_mutation_passes():
    """Testing cases that should pass for enzy_htp.mutation.is_valid_mutation."""
    ref_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    stru = es.PDBParser().get_structure(ref_pdb)
    mutations: List[mut.Mutation] = [
        mut.Mutation(orig='ARG', target='TRP', chain_id='A', res_idx=154),
        mut.Mutation(orig='HIS', target='ALA', chain_id='A', res_idx=201),
        mut.Mutation(orig=None, target='WT', chain_id=None, res_idx=None),
    ]

    for mm in mutations:
        mm.is_valid_mutation(stru)


def test_is_valid_mutation_fails():
    """Testing cases that should fail for enzy_htp.mutation.is_valid_mutation."""
    ref_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    stru = es.PDBParser().get_structure(ref_pdb)
    mutations: List[mut.Mutation] = {
        "wrong data type": [mut.Mutation(orig=1, target='ASP', chain_id='A', res_idx=0), 
                            mut.Mutation(orig="ALA", target='ASP', chain_id='A', res_idx="0"),],
        "empty chain_id": [mut.Mutation(orig='ARG', target='ALA', chain_id='', res_idx=1)],
        "does not exist in structure": [mut.Mutation(orig='ARG', target='ASP', chain_id='B', res_idx=154),
                                        mut.Mutation(orig='ARG', target='ASP', chain_id='A', res_idx=300),],
        "original residue does not match": [mut.Mutation(orig='ALA', target='ASP', chain_id='A', res_idx=154),],
        "unsupported target residue": [mut.Mutation(orig='ARG', target='SEC', chain_id='A', res_idx=154),
                                      mut.Mutation(orig='ARG', target='X', chain_id='A', res_idx=154)],
        "equivalent mutation detected": [mut.Mutation(orig='ARG', target='ARG', chain_id='A', res_idx=154)],
    }

    for msg_finger_p, mm in mutations.items():
        for m in mm:
            with pytest.raises(Exception) as exe:
                m.is_valid_mutation(stru)
            assert exe.type == InvalidMutation
            assert msg_finger_p in exe.value.args[0]


def test_check_repeat_mutation():
    """test to make sure function works as expected"""
    test_mutant = [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "ILE", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]
    assert mut.check_repeat_mutation(test_mutant)
    test_mutant = [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]
    assert not mut.check_repeat_mutation(test_mutant)

def test_remove_repeat_mutation():
    """test to make sure function works as expected"""
    test_mutant = [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "ILE", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]
    assert mut.remove_repeat_mutation(test_mutant, "last") == [
        mut.Mutation("ARG", "ILE", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]
    assert mut.remove_repeat_mutation(test_mutant, "first") == [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]


def test_hash_and_set():
    """test is if __hash__ function is set correct that set() functions correctly"""
    test_mutant = [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]

    test_mutant_2 = [
        mut.Mutation("ARG", "TRP", "A", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("ARG", "TRP", "B", 154),
        mut.Mutation("LYS", "GLN", "A", 37),
    ]
    assert len(set(test_mutant)) == 3
    assert test_mutant[0] == test_mutant[1] # same value
    assert test_mutant[0] is not test_mutant[1] # but different object
    assert set(test_mutant) == set(test_mutant_2)

def test_changed_clone():
    """test to make sure function works as expected"""
    test_mutation = mut.Mutation("ARG", "TRP", "A", 154)
    assert test_mutation.changed_clone(target= "LEU", res_idx= 10) == ("ARG", "LEU", "A", 10)


# == TODO ==
def test_generate_all_mutations():
    """Testing that all possible mutations work for a simple, 1-residue structure."""
    ONE_RES: str = f"{DATA_DIR}/one_res.pdb"
    structure: es.Structure = es.PDBParser().get_structure(ONE_RES)
    all_muts = mut.generate_all_mutations(structure)

    assert len(all_muts[('A', 1)]) == 20

    for mm in all_muts[('A', 1)]:
        assert mm.target != 'P'


def test_size_increase_true():
    """Testing cases where size_increase() should evalue to 'True'"""
    assert mut.size_increase(mut.Mutation(orig='G', target='A', chain_id='',
                                          res_idx=None))
    assert mut.size_increase(mut.Mutation(orig='G', target='S', chain_id='',
                                          res_idx=None))
    assert mut.size_increase(mut.Mutation(orig='ASP', target='P', chain_id='',
                                          res_idx=None))
    assert mut.size_increase(mut.Mutation(orig='Y', target='W', chain_id='',
                                          res_idx=None))


def test_size_increase_false():
    """Testing cases where size_increase() should evalue to 'False'"""
    assert not mut.size_increase(
        mut.Mutation(target='G', orig='A', chain_id='', res_idx=None))
    assert not mut.size_increase(
        mut.Mutation(target='G', orig='S', chain_id='', res_idx=None))
    assert not mut.size_increase(
        mut.Mutation(target='ASP', orig='P', chain_id='', res_idx=None))
    assert not mut.size_increase(
        mut.Mutation(target='Y', orig='W', chain_id='', res_idx=None))
    assert not mut.size_increase(
        mut.Mutation(target='Y', orig='Y', chain_id='', res_idx=None))


def test_size_decrease_true():
    """Testing cases where size_decrease() should evaluate to 'True'"""
    assert mut.size_decrease(mut.Mutation(target='G', orig='A', chain_id='',
                                          res_idx=None))
    assert mut.size_decrease(mut.Mutation(target='G', orig='S', chain_id='',
                                          res_idx=None))
    assert mut.size_decrease(mut.Mutation(target='ASP', orig='P', chain_id='',
                                          res_idx=None))
    assert mut.size_decrease(mut.Mutation(target='Y', orig='W', chain_id='',
                                          res_idx=None))


def test_size_decrease_false():
    """Testing cases where size_decrease() should evaluate to 'False'"""
    assert not mut.size_decrease(
        mut.Mutation(orig='G', target='A', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='G', target='S', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='ASP', target='P', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='Y', target='W', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='Y', target='Y', chain_id='', res_idx=None))


def test_polarity_change_true():
    """Testing cases where polarity_change() should evaluate to 'True'"""
    # negative to positive
    assert mut.polarity_change(
        mut.Mutation(orig='ASP', target='ARG', chain_id='', res_idx=None))
    # positive to negative
    assert mut.polarity_change(
        mut.Mutation(orig='ARG', target='ASP', chain_id='', res_idx=None))
    # positive to neutral
    assert mut.polarity_change(
        mut.Mutation(orig='ARG', target='S', chain_id='', res_idx=None))
    # netural to postiive
    assert mut.polarity_change(
        mut.Mutation(orig='S', target='ARG', chain_id='', res_idx=None))
    # negative to neutral
    assert mut.polarity_change(
        mut.Mutation(orig='ASP', target='S', chain_id='', res_idx=None))
    # neutral to negative
    assert mut.polarity_change(
        mut.Mutation(orig='S', target='ASP', chain_id='', res_idx=None))


def test_polarity_change_false():
    """Testing cases where polarity_change() should evaluate to 'False'"""
    # negative to negative
    assert not mut.polarity_change(
        mut.Mutation(orig='E', target='ASP', chain_id='', res_idx=None))
    # neutral to neutral
    assert not mut.polarity_change(
        mut.Mutation(orig='S', target='T', chain_id='', res_idx=None))
    # positive to positive
    assert not mut.polarity_change(
        mut.Mutation(orig='H', target='ARG', chain_id='', res_idx=None))


def test_same_polarity_true():
    """Testing cases where same_polarity() should evaluate to 'True'"""
    # negative to negative
    assert mut.same_polarity(mut.Mutation(orig='E', target='ASP', chain_id='',
                                          res_idx=None))
    # neutral to neutral
    assert mut.same_polarity(mut.Mutation(orig='S', target='T', chain_id='',
                                          res_idx=None))
    # positive to positive
    assert mut.same_polarity(mut.Mutation(orig='H', target='ARG', chain_id='',
                                          res_idx=None))


def test_same_polarity_false():
    """Testing cases where same_polarity() should evaluate to 'False'"""
    # negative to positive
    assert not mut.same_polarity(
        mut.Mutation(orig='ASP', target='ARG', chain_id='', res_idx=None))
    # positive to negative
    assert not mut.same_polarity(
        mut.Mutation(orig='ARG', target='ASP', chain_id='', res_idx=None))
    # positive to neutral
    assert not mut.same_polarity(
        mut.Mutation(orig='ARG', target='S', chain_id='', res_idx=None))
    # netural to postiive
    assert not mut.same_polarity(
        mut.Mutation(orig='S', target='ARG', chain_id='', res_idx=None))
    # negative to neutral
    assert not mut.same_polarity(
        mut.Mutation(orig='ASP', target='S', chain_id='', res_idx=None))
    # neutral to negative
    assert not mut.same_polarity(
        mut.Mutation(orig='S', target='ASP', chain_id='', res_idx=None))

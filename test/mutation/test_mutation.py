"""Testing the Mutation() namedtuple and associated compatability functions in the enzy_htp.mutation.mutation
submodule.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
        QZ Shao <shaoqz@icloud.com>
Date: 2022-06-15
"""
import logging
import os
from typing import List

from enzy_htp.core.logger import _LOGGER
import enzy_htp.structure as es
import enzy_htp.mutation as mut

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"


def test_decode_mutation_flag_wt(caplog):
    """test if the function behave as expected on WT"""
    assert mut.decode_mutation_flag("WT") == mut.Mutation(
        orig=None, target='WT', chain_id=None, res_idx=None
        )

    original_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    assert mut.decode_mutation_flag("RA154R") == mut.Mutation(
        orig=None, target='WT', chain_id=None, res_idx=None
        )
    assert "equivalent mutation detected" in caplog.text

    _LOGGER.setLevel(original_level)

def test_decode_mutation_flag_default_chainid(caplog):
    """test if the function behave as expected on default chain id"""
    original_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    assert mut.decode_mutation_flag("R154A") == mut.Mutation(
        orig='R', target='A', chain_id='A', res_idx=154
        )
    assert " Using A as default." in caplog.text

    _LOGGER.setLevel(original_level)

def test_is_valid_mutation_passes():
    """Testing cases that should pass for enzy_htp.mutation.is_valid_mutation."""
    mutations: List[mut.Mutation] = [
        mut.Mutation(orig='R', target='D', chain_id='A', res_idx=1),
        mut.Mutation(orig='X', target='D', chain_id='A', res_idx=1),
        mut.Mutation(orig='R', target='D', chain_id='', res_idx=1000),
        mut.Mutation(orig='R', target='D', chain_id=' ', res_idx=1000),
        mut.Mutation(orig='U', target='W', chain_id='A', res_idx=1),
    ]

    for mm in mutations:
        assert mut.is_valid_mutation(mm)


def test_is_valid_mutation_fails():
    """Testing cases that should fail for enzy_htp.mutation.is_valid_mutation."""
    mutations: List[mut.Mutation] = [
        mut.Mutation(orig='R', target='R', chain_id='A', res_idx=1),
        mut.Mutation(orig='R', target='D', chain_id='A', res_idx=0),
        mut.Mutation(orig='R', target='D', chain_id='A', res_idx=0.0),
        mut.Mutation(orig='R', target='D', chain_id=3, res_idx=1000),
        mut.Mutation(orig='R', target='D', chain_id='3', res_idx=1000),
        mut.Mutation(orig='R', target='Z', chain_id='A', res_idx=1),
        mut.Mutation(orig='Z', target='R', chain_id='A', res_idx=1),
        mut.Mutation(orig=1, target='D', chain_id='A', res_idx=0),
        mut.Mutation(orig='R', target=1, chain_id='A', res_idx=0),
        mut.Mutation(orig='R', target='D', chain_id=1, res_idx=0),
        mut.Mutation(orig='R', target='D', chain_id='A', res_idx=False),
    ]

    for mm in mutations:
        assert not mut.is_valid_mutation(mm)


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
    assert mut.size_increase(mut.Mutation(orig='D', target='P', chain_id='',
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
        mut.Mutation(target='D', orig='P', chain_id='', res_idx=None))
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
    assert mut.size_decrease(mut.Mutation(target='D', orig='P', chain_id='',
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
        mut.Mutation(orig='D', target='P', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='Y', target='W', chain_id='', res_idx=None))
    assert not mut.size_decrease(
        mut.Mutation(orig='Y', target='Y', chain_id='', res_idx=None))


def test_polarity_change_true():
    """Testing cases where polarity_change() should evaluate to 'True'"""
    # negative to positive
    assert mut.polarity_change(
        mut.Mutation(orig='D', target='R', chain_id='', res_idx=None))
    # positive to negative
    assert mut.polarity_change(
        mut.Mutation(orig='R', target='D', chain_id='', res_idx=None))
    # positive to neutral
    assert mut.polarity_change(
        mut.Mutation(orig='R', target='S', chain_id='', res_idx=None))
    # netural to postiive
    assert mut.polarity_change(
        mut.Mutation(orig='S', target='R', chain_id='', res_idx=None))
    # negative to neutral
    assert mut.polarity_change(
        mut.Mutation(orig='D', target='S', chain_id='', res_idx=None))
    # neutral to negative
    assert mut.polarity_change(
        mut.Mutation(orig='S', target='D', chain_id='', res_idx=None))


def test_polarity_change_false():
    """Testing cases where polarity_change() should evaluate to 'False'"""
    # negative to negative
    assert not mut.polarity_change(
        mut.Mutation(orig='E', target='D', chain_id='', res_idx=None))
    # neutral to neutral
    assert not mut.polarity_change(
        mut.Mutation(orig='S', target='T', chain_id='', res_idx=None))
    # positive to positive
    assert not mut.polarity_change(
        mut.Mutation(orig='H', target='R', chain_id='', res_idx=None))


def test_same_polarity_true():
    """Testing cases where same_polarity() should evaluate to 'True'"""
    # negative to negative
    assert mut.same_polarity(mut.Mutation(orig='E', target='D', chain_id='',
                                          res_idx=None))
    # neutral to neutral
    assert mut.same_polarity(mut.Mutation(orig='S', target='T', chain_id='',
                                          res_idx=None))
    # positive to positive
    assert mut.same_polarity(mut.Mutation(orig='H', target='R', chain_id='',
                                          res_idx=None))


def test_same_polarity_false():
    """Testing cases where same_polarity() should evaluate to 'False'"""
    # negative to positive
    assert not mut.same_polarity(
        mut.Mutation(orig='D', target='R', chain_id='', res_idx=None))
    # positive to negative
    assert not mut.same_polarity(
        mut.Mutation(orig='R', target='D', chain_id='', res_idx=None))
    # positive to neutral
    assert not mut.same_polarity(
        mut.Mutation(orig='R', target='S', chain_id='', res_idx=None))
    # netural to postiive
    assert not mut.same_polarity(
        mut.Mutation(orig='S', target='R', chain_id='', res_idx=None))
    # negative to neutral
    assert not mut.same_polarity(
        mut.Mutation(orig='D', target='S', chain_id='', res_idx=None))
    # neutral to negative
    assert not mut.same_polarity(
        mut.Mutation(orig='S', target='D', chain_id='', res_idx=None))

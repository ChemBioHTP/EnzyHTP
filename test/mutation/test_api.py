"""Testing the enzy_htp.mutation.api.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
        Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-01-23
"""

import os
from typing import List

import pytest

from enzy_htp.core import file_system as fs
from enzy_htp.core.exception import UnsupportedMethod
from enzy_htp.structure import PDBParser
import enzy_htp.mutation.api as mapi
from enzy_htp.mutation.mutation import Mutation

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
sp = PDBParser()


def test_assign_mutant():
    """test function works as expected using KE07"""
    test_mutation_pattern = (
        "KA162A, {RA154W, HA201A},"
        " {L10A, r:2[resi 254 around 3:all not self]*5}")
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)

    mutants = mapi.assign_mutant(test_stru, test_mutation_pattern)
    # print(mutants)


def test_assign_mutant_dimer():
    """test function works as expected using PuO
    TODO(qz): fix this by filtering polypeptide only"""
    test_mutation_pattern = "r:2[resi 904 around 3:all not self]*5}"
    test_pdb = f"{DATA_DIR}puo_put.pdb"
    test_stru = sp.get_structure(test_pdb)

    mutants = mapi.assign_mutant(test_stru, test_mutation_pattern, chain_sync_list=[("A", "B")], chain_index_mapper={"A": 0, "B": 450})
    print(mutants)


def test_sync_mutation_over_chains():
    """test function works as expected"""
    test_mutants = [[
        Mutation(orig='ARG', target='ALA', chain_id='A', res_idx=3),
        Mutation(orig='ARG', target='TRP', chain_id='A', res_idx=4)
    ], [Mutation(orig='TRP', target='GLY', chain_id='C', res_idx=1),
        Mutation(orig='TRP', target='HIS', chain_id='C', res_idx=2)]]
    test_chain_sync_list = [("A", "B"), ("C", "D")]
    test_chain_index_mapper = {"A": 0, "B": 10, "C": 20, "D": 100}
    result = mapi.sync_mutation_over_chains(test_mutants, test_chain_sync_list, test_chain_index_mapper)
    assert len(result) == 2
    assert len(result[0]) == 4
    assert Mutation(orig='ARG', target='TRP', chain_id='B', res_idx=14) in result[0]
    assert len(result[1]) == 4
    assert Mutation(orig='TRP', target='HIS', chain_id='D', res_idx=82) in result[1]


def test_mutate_stru_with_tleap():
    """test function works as expected"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    mutant = mapi.assign_mutant(test_stru, "R154W")[0]
    mutant_stru = mapi.mutate_stru_with_tleap(test_stru, mutant)

    for new_res, old_res in zip(mutant_stru.residues, test_stru.residues):
        if new_res.idx == 154:
            assert new_res.name == "TRP"
            assert len(new_res.atoms) == 24
        else:
            for new_atom, old_atom in zip(new_res.atoms, old_res.atoms):
                assert new_atom.coord == old_atom.coord
                assert new_atom.name == old_atom.name


def test_check_mutant_stru():
    """test function works as expected"""
    mutant_pdb = f"{DATA_DIR}puo_put_EA323Y_EB773Y_GA171N_GB621N.pdb"
    mutant = [
        Mutation("GLU", "TYR", "A", 323),
        Mutation("GLU", "TYR", "B", 773),
        Mutation("GLY", "ASN", "A", 171),
        Mutation("GLY", "ASN", "B", 621),
    ]
    mutant_stru = sp.get_structure(mutant_pdb)

    mapi.check_mutant_stru(mutant_stru, mutant)


def test_mutate_stru_with_pymol():
    """ test function works as expected """
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    mutant = mapi.assign_mutant(test_stru, "R154W")[0]
    mutant_stru = mapi.mutate_stru_with_pymol(test_stru, mutant)

    for new_res, old_res in zip(mutant_stru.residues, test_stru.residues):
        if new_res.idx == 154:
            assert new_res.name == "TRP"
            assert len(new_res.atoms) == 24
        else:
            for new_atom, old_atom in zip(new_res.atoms, old_res.atoms):
                assert new_atom.name == old_atom.name


def test_mutate_stru_with_pymol_multiple_mutants():
    # test with multiple mutants
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    mutant_2 = mapi.assign_mutant(test_stru, "{R154W, HA201A}")[0]
    mutant_stru_2 = mapi.mutate_stru_with_pymol(test_stru, mutant_2)

    for new_res, old_res in zip(mutant_stru_2.residues, test_stru.residues):
        if new_res.idx == 154:
            assert new_res.name == "TRP"
            assert len(new_res.atoms) == 24
        elif new_res.idx == 201:
            assert new_res.name == "ALA"
        else:
            for new_atom, old_atom in zip(new_res.atoms, old_res.atoms):
                assert new_atom.name == old_atom.name


def test_mutate_stru_with_pymol_multiple_chains():
    # test with > 2 chain enzyme
    test_pdb_2 = f"{DATA_DIR}puo_put.pdb"
    test_stru_2 = sp.get_structure(test_pdb_2)
    mutant_3 = mapi.assign_mutant(test_stru_2, "RB533W")[0]
    mutant_stru_3 = mapi.mutate_stru_with_pymol(test_stru_2, mutant_3)

    for new_res, old_res in zip(mutant_stru_3.residues, test_stru_2.residues):
        if new_res.idx == 533:
            assert new_res.name == "TRP"
        else:
            for new_atom, old_atom in zip(new_res.atoms, old_res.atoms):
                assert new_atom.name == old_atom.name

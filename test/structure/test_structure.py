"""Testing the Structure() class. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
import os
import pytest
import numpy as np
from copy import deepcopy
import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs
from enzy_htp.structure import (
    Structure,
    structure_from_pdb,
    Residue,
    Chain,
    Atom,
    compare_structures,
    merge_right,
)

#TODO(CJ): add tests for Structure.build_ligands()
CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/data/"


def equiv_files(fname1: str, fname2: str, width: int = None) -> bool:
    """Helper method to check if two files are exactly equivalent."""
    def find_first_diff( l1, l2 ):
        for idx, (ch1, ch2) in enumerate(zip(l1,l2)):
            if ch1 != ch2:
                return idx+1
        return -1

    for line_idx, (l1, l2) in enumerate(zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2))):
        if width:
            l1 = l1[:width]
            l2 = l2[:width]
            if len(l2) < width:
                l2 += ' '*(width-len(l2))
            if len(l1) < width:
                l1 += ' '*(width-len(l1))

        first_diff = find_first_diff(l1,l2)
        if first_diff != -1:
            diff = [' ']*min(len(l1),len(l2))
            diff[first_diff] = '^'
            
            print(f"'{l1}'")
            print(f"'{l2}'")
            print(''.join(diff))
            print(f"Difference encountered on line {line_idx}. '{fname1}' and '{fname2}' are NOT equivalent")
            return False
    return True


def test_load_structure():
    """Testing that the structure_from_pdb() method works for a basic example."""
    raw_answer = [
        ("A", "THR", 1),
        ("A", "THR", 2),
        ("A", "CYS", 3),
        ("A", "CYS", 4),
        ("A", "PRO", 5),
        ("A", "SER", 6),
        ("A", "ILE", 7),
        ("A", "VAL", 8),
        ("A", "ALA", 9),
        ("A", "ARG", 10),
        ("A", "SER", 11),
        ("A", "ASN", 12),
        ("A", "PHE", 13),
        ("A", "ASN", 14),
        ("A", "VAL", 15),
        ("A", "CYS", 16),
        ("A", "ARG", 17),
        ("A", "LEU", 18),
        ("A", "PRO", 19),
        ("A", "GLY", 20),
        ("A", "THR", 21),
        ("A", "PRO", 22),
        ("A", "SER", 22),
        ("A", "GLU", 23),
        ("A", "ALA", 24),
        ("A", "LEU", 25),
        ("A", "ILE", 25),
        ("A", "CYS", 26),
        ("A", "ALA", 27),
        ("A", "THR", 28),
        ("A", "TYR", 29),
        ("A", "THR", 30),
        ("A", "GLY", 31),
        ("A", "CYS", 32),
        ("A", "ILE", 33),
        ("A", "ILE", 34),
        ("A", "ILE", 35),
        ("A", "PRO", 36),
        ("A", "GLY", 37),
        ("A", "ALA", 38),
        ("A", "THR", 39),
        ("A", "CYS", 40),
        ("A", "PRO", 41),
        ("A", "GLY", 42),
        ("A", "ASP", 43),
        ("A", "TYR", 44),
        ("A", "ALA", 45),
        ("A", "ASN", 46),
    ]
    fixed_answer = list(
        map(lambda pr: (pr[0], chem.convert_to_one_letter(pr[1]), pr[2]), raw_answer)
    )
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    struct: Structure = structure_from_pdb(TEST_FILE)
    assert struct
    assert struct.residue_state == fixed_answer


def test_structure_ctor_bad_input():
    """Testing that the Structure() ctor fails when given duplicate chains."""
    chain1 = Chain("A", [Residue('A.XXX.1', [Atom(x_coord=1, y_coord=2, z_coord=3)])])
    chain2 = Chain("B", [Residue('B.XXX.1', [Atom(x_coord=2, y_coord=2, z_coord=3)])])
    chain3 = Chain("A", [Residue('A.XXX.1', [Atom(x_coord=1, y_coord=2, z_coord=3)])])

    with pytest.raises(SystemExit) as exe:
        struct = Structure([chain1, chain2, chain3])

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_residue_state():
    """Checking that the Structureresidue_state method properly returns the residue state for the Structure()."""
    start = [("A", "ARG", 1), ("A", "HIS", 2), ("A", "ARG", 3), ("A", "HIS", 4)]
    answer = [("A", "R", 1), ("A", "H", 2), ("A", "R", 3), ("A", "H", 4)]
    keys = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start))
    chain1 = Chain("A", list(map(lambda kk: Residue(kk, list()), keys)))
    ss = Structure(
        [chain1]
    )  # TODO(CJ): need error in Structure() ctor if NOT a list of chains
    assert ss.residue_state == answer


def test_structure_same_sequence():
    """Ensuring that the Chain.is_same_sequence() method works."""
    start1 = [("A", "ARG", 1), ("A", "HIS", 2), ("A", "ARG", 3), ("A", "HIS", 4)]
    start2 = [("B", "ARG", 1), ("B", "HIS", 2), ("B", "ARG", 3), ("B", "HIS", 4)]
    start3 = [("A", "ARG", 1), ("A", "ARG", 2), ("A", "ARG", 3), ("A", "HIS", 4)]
    start4 = [("A", "ARG", 1), ("A", "HIS", 2), ("A", "ARG", 4), ("A", "HIS", 4)]
    keys1 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start1))
    keys2 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start2))
    keys3 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start3))
    keys4 = list(map(lambda pr: f"{pr[0]}.{pr[1]}.{pr[2]}", start4))
    # testing the same chain
    chain1 = Chain("A", list(map(lambda kk: Residue(kk, list()), keys1)))
    chain2 = Chain("A", list(map(lambda kk: Residue(kk, list()), keys1)))
    assert chain1.is_same_sequence(chain2)
    assert chain2.is_same_sequence(chain1)
    # testing different chain name
    chain2 = Chain("B", list(map(lambda kk: Residue(kk, list()), keys2)))
    assert not chain1.is_same_sequence(chain2)
    assert not chain2.is_same_sequence(chain1)
    # testing different residue name
    chain2 = Chain("A", list(map(lambda kk: Residue(kk, list()), keys3)))
    assert not chain1.is_same_sequence(chain2)
    assert not chain2.is_same_sequence(chain1)
    # testing different residue number
    chain2 = Chain("A", list(map(lambda kk: Residue(kk, list()), keys4)))
    assert not chain1.is_same_sequence(chain2)
    assert not chain2.is_same_sequence(chain1)


def test_compare_structures_equiv():
    """Ensuring the enzy_htp.structure.compare_structures() free function identifies two identical structures."""
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE)
    structure2: Structure = structure_from_pdb(TEST_FILE)
    result = compare_structures(structure1, structure2)
    assert result == {"left": [], "right": []}


def test_compare_structures_not_equiv():
    """Ensuring the enzy_htp.structure.compare_structures() free functions identifies differences between functions."""
    # TODO(CJ): include non-canonical residues in here
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE)
    structure2: Structure = structure_from_pdb(TEST_FILE)
    res_to_remove1: Residue = structure1.chains[0][-1]
    res_to_remove2: Residue = structure2.chains[0][0]
    # one missing residue
    del structure1.chains[0][-1]
    result = compare_structures(structure1, structure2)
    assert result == {"left": [], "right": [res_to_remove1.residue_key]}
    # two missing residues
    del structure2.chains[0][0]
    result = compare_structures(structure1, structure2)
    assert result == {
        "left": [res_to_remove2.residue_key],
        "right": [res_to_remove1.residue_key],
    }


def test_merge_right_canonical_only(): #@shaoqz: @imp a better idea is to merge 2 different protein.
    """Merges differences between two Structure() objects containing only canonical residues."""
    # TODO(CJ): include non-canonical residues in here
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE)
    structure2: Structure = structure_from_pdb(TEST_FILE)
    res_to_remove1: Residue = structure2.chains[0][-1]
    del structure2.chains[0][-1]
    result = compare_structures(structure1, structure2)
    assert result == {"left": [res_to_remove1.residue_key], "right": []}
    result = merge_right(structure1, structure2)
    assert structure1 == result


def test_merge_right_with_ligand():
    """Merges differences between two Structure() objects containing a Ligand() object."""
    TEST_FILE = f"{DATA_DIR}/FAcD-FA-ASP_rmW.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE)
    structure2: Structure = structure_from_pdb(TEST_FILE)
    structure2.remove_chain("B")
    assert structure1 != structure2
    structure2: Structure = merge_right(structure1, structure2)
    assert structure1 == structure2
    print(DATA_DIR)
    # assert False


def test_merge_right_with_metal_atom():
    """Merges differences between two Structure() objects containing a MetalUnit() object."""
    TEST_FILE = f"{DATA_DIR}/1NVG.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE)
    structure2: Structure = structure_from_pdb(TEST_FILE)
    # NOTE(CJ): Below is NOT the recommended way to manipulate this stuff
    del structure2._chains[0]._residues[-1]
    print(structure1._chains[0]._residues)
    print(structure2._chains[0]._residues)
    assert structure1 != structure2
    structure2: Structure = merge_right(structure1, structure2)
    assert structure1 == structure2


def test_round_trip_pdb():
    """Ensuring that the Structure() class be loaded into a .pdb and saved back in a round trip without error."""
    # FIXME(CJ): This test doesn't currently work for 1NVG: figure out the PDBline stuff
    # that will make this work
    TEST_FILE = f"{DATA_DIR}/1NVG.pdb"
    actual_file = f"{DATA_DIR}/1NVG_cpy.pdb"
    structure1: Structure = structure_from_pdb(TEST_FILE, 'all')
    fs.safe_rm(actual_file)
    assert not os.path.exists(actual_file)
    structure1.to_pdb(actual_file)
    assert os.path.exists(actual_file)
    assert equiv_files(TEST_FILE, actual_file, 60)
    fs.safe_rm(actual_file)
    assert not os.path.exists(actual_file)


def test_atoms(): # TODO(shaoqz) wait for test
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    struct: Structure = structure_from_pdb(TEST_FILE)
    print(struct.atoms)

"""Testing the Structure() class. 
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2022-04-03
"""
import os
import time
import pytest
import numpy as np
from copy import deepcopy
import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.structure import (
    Structure,
    Residue,
    Chain,
    Atom,
    Ligand,
)

#TODO(CJ): add tests for Structure.build_ligands()
CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/data/"
sp = PDBParser()

def equiv_files(fname1: str, fname2: str, width: int = None) -> bool:
    """Helper method to check if two files are exactly equivalent."""

    def find_first_diff(l1, l2):
        for idx, (ch1, ch2) in enumerate(zip(l1, l2)):
            if ch1 != ch2:
                return idx + 1
        return -1

    for line_idx, (l1, l2) in enumerate(zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2))):
        if width:
            l1 = l1[:width]
            l2 = l2[:width]
            if len(l2) < width:
                l2 += " " * (width - len(l2))
            if len(l1) < width:
                l1 += " " * (width - len(l1))

        first_diff = find_first_diff(l1, l2)
        if first_diff != -1:
            diff = [" "] * min(len(l1), len(l2))
            diff[first_diff] = "^"

            print(f"'{l1}'")
            print(f"'{l2}'")
            print("".join(diff))
            print(f"Difference encountered on line {line_idx}. '{fname1}' and '{fname2}' are NOT equivalent")
            return False
    return True


def test_structure_ctor_bad_input():
    """Testing that the Structure() ctor fails when given duplicate chains."""
    Atom_1 = Atom.from_biopandas({"atom_name": "X", "x_coord": 1, "y_coord": 2, "z_coord": 3})
    Atom_2 = Atom.from_biopandas({"atom_name": "X", "x_coord": 3, "y_coord": 2, "z_coord": 1})
    chain1 = Chain("A", [Residue(1, "XXX", [Atom_1])])
    chain2 = Chain("B", [Residue(2, "XXX", [Atom_2])])
    chain3 = Chain("A", [Residue(3, "XXX", [Atom_1])])

    with pytest.raises(SystemExit) as exe:
        struct = Structure([chain1, chain2, chain3])

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_residues():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    assert len(stru.residues) == 580


def test_ligands():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.ligands) == 2
    assert stru.ligands[0].name == "4CO"


def test_solvents():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.solvents) == 296
    assert stru.solvents[0].name == "HOH"


def test_peptides():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.polypeptides) == 2
    assert tuple(map(lambda x: x.name, stru.polypeptides)) == ("A", "B")


def test_sequence():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    target_seq: str = "ATGGNLPDVASHYPVAYEQTLDGTVGFVIDEMTPERATASVEVTDTLRQRWGLVHGGAYCALAEMLATEATVAVVHEKGMMAVGQSNHTSFFRPVKEGHVRAEAVRIHAGSTTWFWDVSLRDDAGRLCAVSSMSIAVRPRRD"
    assert stru.sequence["A"] == target_seq


def test_is_idx_subset_subset():
    """test the checker for idx subset with a handmade subset stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(subset_pdb_file_path)
    assert stru.is_idx_subset(target_stru)


def test_is_idx_subset_non_res_subset():
    """test the checker for idx subset with a handmade non-subset (residue) stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    non_subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_non_seq_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(non_subset_pdb_file_path)
    assert not stru.is_idx_subset(target_stru)


def test_is_idx_subset_non_ch_subset():
    """test the checker for idx subset with a handmade non-subset (chain id) stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    non_subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_non_chain_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(non_subset_pdb_file_path)
    assert not stru.is_idx_subset(target_stru)


def test_atoms():  # TODO(shaoqz) wait for test
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    struct: Structure = PDBParser().get_structure(TEST_FILE)
    assert struct.atoms


def test_deepcopy():
    """test the hehavior of copy.deepcopy on Structure()
    context"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    new_stru = deepcopy(stru)
    # ensure the list is new
    assert id(stru) != id(new_stru)
    # ensure children are pointing to the parent
    for ch in new_stru.chains:
        assert ch.parent is new_stru
        for res in ch:
            assert res.parent is ch
            for atom in res:
                assert atom.parent is res


def test_deepcopy_profiling():
    """Test the performance and correctness of deepcopy on Structure objects.
    
    This test evaluates the time cost of deepcopy operation on a large solvated structure
    and verifies that the deep copy maintains correct parent-child relationships.
    """
    # Load a relatively large structure for performance testing
    stru = PDBParser().get_structure(f"{DATA_DIR}test_pdb_parser_solvated.pdb")
    
    # Record structure statistics before copying
    num_atoms = len(stru.atoms)
    
    # Perform deepcopy with timing
    start_time = time.perf_counter()
    new_stru = deepcopy(stru)
    end_time = time.perf_counter()
    
    deepcopy_time = end_time - start_time
    print(deepcopy_time)

    max_time_per_1000_atoms = 0.2  # 200ms per 1000 atoms (generous threshold)
    expected_max_time = (num_atoms / 1000) * max_time_per_1000_atoms
    assert deepcopy_time < expected_max_time, f"Deepcopy took too long: {deepcopy_time:.4f}s > {expected_max_time:.4f}s (expected < {expected_max_time:.4f}s)"


def test_clone():
    """test the behavior of clone on Structure() context"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    new_stru = stru.clone()
    # ensure the list is new
    assert id(stru) != id(new_stru)
    # ensure children are pointing to the parent
    for ch in new_stru.chains:
        assert ch.parent is new_stru
        for res in ch:
            assert res.parent is ch
            for atom in res:
                assert atom.parent is res
    # ensure there are no shared objects but the same topology
    for och, nch in zip(stru.chains, new_stru.chains):
        assert och.name == nch.name
        assert id(och) != id(nch)
        for ores, nres in zip(och.residues, nch.residues):
            assert ores.idx == nres.idx
            assert id(ores) != id(nres)
            for oatom, natom in zip(ores.atoms, nres.atoms):
                assert oatom.name == natom.name
                assert id(oatom) != id(natom)


def test_clone_profile():
    """test the behavior of clone on Structure() context"""
    # Load a relatively large structure for performance testing
    stru = PDBParser().get_structure(f"{DATA_DIR}test_pdb_parser_solvated.pdb")
    
    # Record structure statistics before copying
    num_atoms = len(stru.atoms)
    
    # Perform deepcopy with timing
    start_time = time.perf_counter()
    new_stru = stru.clone()
    end_time = time.perf_counter()
    
    deepcopy_time = end_time - start_time
    print(deepcopy_time)

    max_time_per_1000_atoms = 0.2  # 200ms per 1000 atoms (generous threshold)
    expected_max_time = (num_atoms / 1000) * max_time_per_1000_atoms
    assert deepcopy_time < expected_max_time, f"Deepcopy took too long: {deepcopy_time:.4f}s > {expected_max_time:.4f}s (expected < {expected_max_time:.4f}s)"


def test_find_residue_with_key(caplog): # TODO caplog does not work when logger is not propagate. fix them by a context manager.
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    res_a47 = stru.find_residue_with_key(("A", 47))
    assert res_a47.name == "THR"
    # warning case
    res_a47.idx = 48
    stru.find_residue_with_key(("A", 48))
    assert "More than 1 residue with key: ('A', 48)" in caplog.text
    # non case
    assert not stru.find_residue_with_key(("A", 448))
    assert "Didn't find any residue with key: ('A', 448)" in caplog.text


def test_chemical_diversity_ligand():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert stru.chemical_diversity == {'polypeptide': 2, 'ligand': {'N', 'S', 'O', 'P', 'C'}, 'solvents': {'O'}}

def test_chemical_diversity_metal():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1NVG.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert stru.chemical_diversity == {'polypeptide': 1, 'metalcenters': {'Zn'}}

def test_chemical_diversity_nucleic_acid():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}cas9_4oo8.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert stru.chemical_diversity == {'polypeptide': 2, 'ligand': {'N', 'C', 'O', 'P'}, 'solvents': {'O'}, 'nucleic_acid': ''}

def test_add_chain():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    new_chain = Chain("X", [], None)
    stru.add(new_chain)
    assert stru.chain_mapper.get("X", False) is not False

def test_add_residue():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    new_res = Ligand(
        residue_idx=0,
        residue_name="LIG",
        atoms=[],
        parent=None,)
    stru.add(new_res)
    assert stru["E"].residues[0] is new_res

def test_get(caplog):
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert stru.get("A") is stru["A"]
    assert stru.get("A.10") is stru["A"].find_residue_idx(10)
    assert stru.get("A.10.CA") is stru["A"].find_residue_idx(10).find_atom_name("CA")
    with EnablePropagate(_LOGGER):
        with pytest.raises(ValueError) as e:
            stru.get("A.1")
        assert "Unable to locate " in caplog.text
        with pytest.raises(ValueError) as e:
            stru.get("A.10.X")
        assert "Unable to locate " in caplog.text

def test_residue_index_mapper(): # TODO
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    raise Exception("TODO")

def test_is_same_topology():
    """as name"""
    test_other = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S_geom_1.pdb")
    test_self = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    assert test_self.is_same_topology(test_other)

def test_hydrogens():
    """as nam. result verified by pymol"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    assert len(test_stru.hydrogens()) == 2000

def test_hydrogens_polypeptide_only():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    assert len(test_stru.hydrogens(polypeptide_only=True)) == 1996

def test_net_chrgspin_mapper():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    assert test_stru.ncaa_chrgspin_mapper == {"H5J" : (0,1)}

def test_amino_acids():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    
    assert len(test_stru.amino_acids) == 253


def test_is_topology_subset_atomic_identical():
    """Test is_topology_subset_atomic with identical structures"""
    test_stru = sp.get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    cloned_stru = test_stru.clone(with_connectivity=False)  # Clone without connectivity for clean test
    
    # Identical structure should be a subset of itself
    assert test_stru.is_topology_subset_atomic(test_stru)
    
    # Cloned structure should be a subset of original (same topology)
    assert test_stru.is_topology_subset_atomic(cloned_stru)
    assert cloned_stru.is_topology_subset_atomic(test_stru)


def test_is_topology_subset_atomic_partial():
    """Test is_topology_subset_atomic with partial structures"""
    test_stru = sp.get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    
    # Create a smaller structure with just the first chain
    first_chain = test_stru.chains[0]
    partial_stru = Structure([first_chain.clone(is_clone_root=True)])
    
    # Partial structure should be a subset of full structure
    assert test_stru.is_topology_subset_atomic(partial_stru)
    
    # Full structure should NOT be a subset of partial structure
    assert not partial_stru.is_topology_subset_atomic(test_stru)


def test_create_atom_mapping_identical():
    """Test create_atom_mapping with identical structures"""
    test_stru = sp.get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    cloned_stru = test_stru.clone(with_connectivity=False)
    
    # Create mapping between identical structures
    mapping = test_stru.create_atom_mapping(cloned_stru)
    
    # Should have mapping for every atom
    assert len(mapping) == len(cloned_stru.atoms)
    assert len(mapping) == len(test_stru.atoms)
    
    # Verify mapping is correct by checking atom keys
    for cloned_atom, original_atom in mapping.items():
        assert cloned_atom.key == original_atom.key


def test_create_atom_mapping_partial():
    """Test create_atom_mapping with partial structure"""
    test_stru = sp.get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    
    # Create partial structure with first chain only
    first_chain = test_stru.chains[0]
    partial_stru = Structure([first_chain.clone(is_clone_root=False)])
    
    # Create mapping from partial to full structure
    mapping = test_stru.create_atom_mapping(partial_stru)
    
    # Should have mapping for every atom in partial structure
    assert len(mapping) == len(partial_stru.atoms)
    assert len(mapping) < len(test_stru.atoms)
    
    # Verify mapping correctness
    for partial_atom, original_atom in mapping.items():
        assert partial_atom.key == original_atom.key

"""Testing the connectivity API in enzy_htp.structure.structure_enhancement
Author: Sebastian Stull <sebastian.l.stull@vanderbilt.edu>
Date: 2025-02-13
"""

from collections import defaultdict
import os
import copy
import pytest

from enzy_htp.preparation.clean import remove_solvent
import enzy_htp.structure as struct
from enzy_htp.structure.structure_enchantment import connectivity
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER
import enzy_htp.structure.structure_region as stru_regi
from enzy_htp.structure.structure_region import capping
import enzy_htp.structure.structure_selection as stru_sele


CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
NCAA_LIB = f"{CURR_DIR}/../ncaa_lib"


def test_connectivity_maa():
    test_stru = struct.PDBParser().get_structure(f"{DATA_DIR}/3FCR_modified.pdb")
    test_stru.assign_ncaa_chargespin({"LLP": (-2, 1)})
    remove_solvent(test_stru)

    connectivity.init_connectivity(test_stru)

    fs.safe_rm(f"{NCAA_LIB}/LLP_any.prepin")

    assert test_stru.modified_residue[0].is_connected()


def test_connected_structure_deepcopy():
    """Test to make sure deepcopy works on connected structures.
    NOTE: This test used to fail because the default deep copier walks the entire
    connectivity graph depth-first, so on large molecules it can exceed Python's
    recursion limit and raise RecursionError (even though cycles are memoized).
    This should be fixed now by the change in DoubleLinkedNode on 2025-08-05."""
    # Load the same structure that causes the error in integration tests
    test_stru = struct.PDBParser().get_structure(f"{DATA_DIR}/3FCR_modified.pdb")
    test_stru.assign_ncaa_chargespin({"LLP": (-2, 1), "RLP": (-2, 1)})
    remove_solvent(test_stru)
    
    # Test deepcopy BEFORE connectivity initialization
    copied_stru_before = copy.deepcopy(test_stru)
    
    # Now initialize connectivity (this creates the circular references)
    connectivity.init_connectivity(test_stru)
    
    # Test deepcopy AFTER connectivity initialization
    copied_stru_after = copy.deepcopy(test_stru)
    import pdb;pdb.set_trace()

    # Verify the copied structure maintains connectivity
    # Test that the modified residue is still connected in the copied structure
    assert len(copied_stru_after.modified_residue) > 0, "Copied structure should have modified residue"
    assert copied_stru_after.modified_residue[0].is_connected(), "Copied structure's modified residue should be connected"
    
    # Verify connectivity integrity by checking a few atoms from the modified residue
    # Get the first modified residue (LLP) and check connectivity of its atoms
    original_maa = test_stru.modified_residue[0]
    copied_maa = copied_stru_after.modified_residue[0]
    
    # Test connectivity for 5 key atoms: N, CA, C, CB, and one sidechain atom
    test_atom_names = ['N', 'CA', 'C', 'CB', 'P1']  # P1 is specific to LLP residue
    
    for atom_name in test_atom_names:
        original_atom = None
        copied_atom = None
        
        # Find the atoms in both structures
        for atom in original_maa.atoms:
            if atom.name == atom_name:
                original_atom = atom
                break
        for atom in copied_maa.atoms:
            if atom.name == atom_name:
                copied_atom = atom
                break
        
        # If atom exists in original, it should exist in copy with same connectivity
        if original_atom is not None:
            assert copied_atom is not None, f"Atom {atom_name} should exist in copied structure"
            assert original_atom.is_connected() == copied_atom.is_connected(), \
                f"Atom {atom_name} connectivity should be preserved in copy"
            
            # Check that connected atoms count is the same
            original_connections = len(original_atom.connect_atoms)
            copied_connections = len(copied_atom.connect_atoms)
            assert original_connections == copied_connections, \
                f"Atom {atom_name} should have same number of connections ({original_connections}) in copy, got {copied_connections}"

def test_connected_structure_deepcopy_ref_leak():
    """Test to make sure deepcopy works on connected structures and there are no reference to the original structure in the copied structure."""
    # Load the same structure that causes the error in integration tests
    test_stru = struct.PDBParser().get_structure(f"{DATA_DIR}/3FCR_modified.pdb")
    test_stru.assign_ncaa_chargespin({"LLP": (-2, 1), "RLP": (-2, 1)})
    remove_solvent(test_stru)
    # Now initialize connectivity (this creates the circular references)
    connectivity.init_connectivity(test_stru)
    
    # Test deepcopy AFTER connectivity initialization
    copied_stru = copy.deepcopy(test_stru)

    original_atom_ids = set()

    for atom in test_stru.atoms:
        original_atom_ids.add(id(atom))

    assert original_atom_ids

    copied_atom_ids = set()

    for atom in copied_stru.atoms:
        copied_atom_ids.add(id(atom))

    assert not (original_atom_ids & copied_atom_ids)

    leaked_refs = defaultdict(list)
    for atom in copied_stru.atoms:
        connect_list = atom.connect_atoms
        for nb_atom in connect_list:
            if id(nb_atom) in original_atom_ids:
                leaked_refs[atom.key].append(nb_atom.key)
    leak_msg = (
        f"found {len(leaked_refs)} connections still pointing to original structure atoms.\n"
        "Example:"
        )
    for i, (atom_key, nb_atom_keys) in enumerate(leaked_refs.items()):
        if i >= 5:
            break
        leak_msg += f"\n    {atom_key} -cnt->"
        for nb_atom_key in nb_atom_keys:
            leak_msg += f"\n        {nb_atom_key} (Old)"
    assert not leaked_refs, leak_msg


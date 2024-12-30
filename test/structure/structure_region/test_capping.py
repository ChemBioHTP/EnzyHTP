"""Testing the enzy_htp.structure.structure_region.capping.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-31
"""
import pytest
import os
from copy import deepcopy
from typing import List

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.core.file_system as fs
import enzy_htp.structure.structure_region.capping as capping
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_selection as stru_sele
from enzy_htp.structure.structure_region.residue_caps import SUPPORTED_CAPS, ResidueCap, CH3Cap 
from enzy_htp.structure.structure_io.xyz_io import XYZParser

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
sp = PDBParser()
xyzp = XYZParser()

def test_capping_with_residue_terminals_tri_a(helpers):
    """as name. use an example region from selection.
    Just make sure not execptions is raised for now.
    TODO complete this test"""
    test_stru = sp.get_structure(f"{DATA_DIR}tri_alanine.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2")
    test_region = stru_regi.StructureRegion(atoms=sele.atoms)
    unique_names = set()
    for atom in test_region.atoms:
        unique_names.add(atom.name)

    assert len(test_region.atoms) == len(unique_names)
        

def test_capping_with_residue_terminals_cap_types_supported():
    """Checking that all the cap types that are supposed to work do not throw errors using capping_with_residue_terminals().
    This code basically back checks the methods listed in a dict() so if one is not listed there, this test could pass but 
    be incorrect.."""
    test_stru = sp.get_structure(f"{DATA_DIR}tri_alanine.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2")

    cap_types:List[str] = list(SUPPORTED_CAPS.keys())

    for nterm_cap in cap_types:
        for cterm_cap in cap_types:
            test_region = stru_regi.StructureRegion(atoms=sele.atoms)
            try:                
                capping.capping_with_residue_terminals(test_region, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
            except Exception as exc:
                assert False, f"capping_with_residue_terminals() failed with {nterm_cap=} {cterm_cap=}"




def test_capping_with_residue_terminals_cap_types_unique_names():
    """Checking that all the cap types do not conflict with any other atom names of other caps from the other terminus 
    for capping_with_residue_terminals().
    Atom names can be shared between different cterm caps, but no names can be shared between a cterm and nterm cap
    under ANY circumstances."""
    test_stru = sp.get_structure(f"{DATA_DIR}tri_alanine.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2")

    cap_types:List[str] = list(SUPPORTED_CAPS.keys())

    for nterm_cap in cap_types:
        for cterm_cap in cap_types:
            test_region = stru_regi.StructureRegion(atoms=sele.atoms)
            capping.capping_with_residue_terminals(test_region, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
   
            unique_names = set()
            for atom in test_region.atoms:
                unique_names.add(atom.name)
            assert len(test_region.atoms) == len(unique_names)



def test_capping_with_residue_terminals_realistic_bond_distances():
    """Checking that all the cap types from capping_with_residue_terminals() have realistic bond distances with their link
    Residue()'s and Atom()'s. This is a basic sanity check
    without invoking a bunch of math that ensures the ResidueCap()'s are being put at reasonable distances
    from the original Residue()'s."""
    test_stru = sp.get_structure(f"{DATA_DIR}tri_alanine.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2")

    cap_types:List[str] = list(SUPPORTED_CAPS.keys())

    for nterm_cap in cap_types:
        for cterm_cap in cap_types:
            test_region = stru_regi.StructureRegion(atoms=sele.atoms)
            capping.capping_with_residue_terminals(test_region, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
         
            for atom in test_region.atoms:
                if not isinstance(atom.parent, ResidueCap):
                    continue
                link_distance:float=atom.parent.link_atom.distance_to(atom.parent.socket_atom)
                assert link_distance >= 1.0 and link_distance <= 1.50



def test_capping_with_residue_terminals_correct_num_caps():
    """Testing that the correct number of caps are made in capping_with_residue_terminals() for a variety of setups
     and edge cases, i.e. the inclusion
    of an already n-terminal or c-terminal Residue()."""

    test_stru = sp.get_structure(f"{DATA_DIR}tri_alanine.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 1+2+3")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 0


    sele = stru_sele.select_stru(
        test_stru, "resi 1")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 1 


    sele = stru_sele.select_stru(
        test_stru, "resi 3")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 1 


    sele = stru_sele.select_stru(
        test_stru, "resi 1+2")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 1 

    sele = stru_sele.select_stru(
        test_stru, "resi 2")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 2 


    sele = stru_sele.select_stru(
        test_stru, "resi 1+3")

    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    caps = set()

    capping.capping_with_residue_terminals(test_region)

    for atom in test_region.atoms:
        if not isinstance(atom.parent, ResidueCap):
            continue
        caps.add(atom.parent)

    assert len(caps) == 2 

def test_capping_with_residue_terminals(helpers):
    """as name. use an example region from selection.
    Just make sure not execptions is raised for now.
    TODO complete this test"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 1+2+3+7+8")
    test_region = stru_regi.StructureRegion(atoms=sele.atoms)
    capping.capping_with_residue_terminals(test_region, nterm_cap='CH3', cterm_cap='CH3')
    # test
    atoms = test_region.atoms
    test_file = f"{WORK_DIR}test_capping_1.xyz"
    answer_file = f"{DATA_DIR}answer_capping_1.xyz"
    lines = [
        str(len(atoms)),
        "",]
    for at in atoms:
        lines.append(f"{at.element} {at.coord[0]} {at.coord[1]} {at.coord[2]}")
    fs.write_lines(test_file, lines)
    assert helpers.equiv_files(test_file, answer_file, consider_order=False)

    fs.safe_rm(test_file)

def test_capping_modaa(helpers):
    # caps a modaa with OH and H at the cterm and nterm, respectively
    test_stru = sp.get_structure(f"{DATA_DIR}3FCR_modified.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 288")
    test_region = stru_regi.StructureRegion(atoms=sele.atoms)

    cterm_cap = "OH"
    nterm_cap = "H"

    capping.capping_with_residue_terminals(test_region, nterm_cap=nterm_cap, cterm_cap=cterm_cap)

    answer_file = f"{DATA_DIR}answer_capping_3.xyz"
    
    result = xyzp.get_file_str(test_region)

    assert fs.content_from_file(answer_file) == result


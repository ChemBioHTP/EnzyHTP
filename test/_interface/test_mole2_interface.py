"""Testing enzy_htp._interface.mole2_interface.py
Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Date: 2025-05-18
"""

from os import path
from typing import List, Tuple
from os import path
from typing import List
from enzy_htp import interface, _LOGGER, PDBParser
from enzy_htp._interface import Mole2Cavity
import pyvista as pv

sp = PDBParser()
mole2 = interface.mole2

DATA_DIR = path.join(path.dirname(__file__), "data")
WORK_DIR = path.join(path.dirname(__file__), "work_dir")

def test_identify_cavities():
    """Test the `interface.mole2.identify_cavities` function."""
    pdb_filepath = path.join(DATA_DIR, "cavity_calc", "aclHMT-ETI-SAH_no-ETI.pdb")
    cavities: List[Mole2Cavity] = interface.mole2.identify_cavities(pdb_filepath, work_dir=WORK_DIR)
    
    has_target_cavity = False
    for i, cavity in enumerate(cavities):
        _LOGGER.info(f"Cavity {i+1}, Open Edges: {cavity.mesh_.n_open_edges}, Volume: {cavity.volume()}")
        if (cavity.volume() > 900):
            has_target_cavity = True
        else:
            pass
        continue
    assert has_target_cavity

def test_parse_cavity():
    """Test the `interface.mole2._parse_cavity` function."""
    mesh_filepath = path.join(DATA_DIR, "cavity_calc", "cavity_1.mesh")
    cavity = interface.mole2._parse_cavity(
        mesh_filepath=mesh_filepath, 
        probe=interface.mole2.config_.PROBE, 
        inner=interface.mole2.config_.INNER,
        mesh_density=interface.mole2.config_.MESH_DENSITY,
    )
    # print(cavity.volume())
    assert abs(cavity.volume() - 946.993) < 1   # Inconsistent result comparing with Mole2 output.

def test_read_cavity_from_xml():
    """Test the `interface.mole2._read_cavity_from_xml` function."""
    xml_filepath = path.join(DATA_DIR, "cavity_calc", "cavities.xml")
    cavity_id = 3  # Use the 3rd cavity for test.
    
    volume, boundary_residue_keys, inner_residue_keys = interface.mole2._read_cavity_from_xml(xml_filepath, cavity_id)
    
    # Verify value type.
    assert isinstance(boundary_residue_keys, list)
    assert isinstance(inner_residue_keys, list)
    
    # Verify content.
    expected_boundary_keys = [("A", 55), ("A", 58), ("A", 63), ("A", 64), ("A", 126), ("A", 127), ("A", 130)]
    expected_inner_keys = [("A", 50), ("A", 54), ("A", 61), ("A", 62), ("A", 92), ("A", 94), ("A", 95), ("A", 97), ("A", 125), ("A", 131), ("A", 132), ("A", 160), ("A", 162), ("A", 195)]
    # expected_boundary = "ILE 55 A,ARG 58 A,GLY 63 A,PRO 64 A,ARG 126 A,ASP 127 A,ILE 130 A"
    # expected_inner = "LEU 50 A,LEU 54 A,ILE 61 A,GLY 62 A,LEU 92 A,SER 94 A,PHE 95 A,TYR 97 A,VAL 125 A,GLY 131 A,GLN 132 A,VAL 160 A,TYR 162 A,ILE 195 A"
    
    assert abs(volume - 415.104) < 0.1
    assert boundary_residue_keys == expected_boundary_keys
    assert inner_residue_keys == expected_inner_keys
    
    # Test nonexistent cavity_id.
    try:
        interface.mole2._read_cavity_from_xml(xml_filepath, 999)
        assert False, "Should raise ValueError for non-existent cavity_id"
    except ValueError:
        pass

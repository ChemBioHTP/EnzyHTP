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
        mesh_filepath, 
        probe=interface.mole2.config_.PROBE, 
        inner=interface.mole2.config_.INNER,
        mesh_density=interface.mole2.config_.MESH_DENSITY,
    )
    # print(cavity.volume())
    assert abs(cavity.volume() - 946.993) < 1   # Inconsistent result comparing with Mole2 output.
    
    # cavity.mesh_.plot(
    #     # scalars=np.array(3),
    #     # cpos=[-1, 1, 0.5],
    #     show_scalar_bar=False,
    #     show_edges=True,
    #     line_width=5,
    # )
    # plotter = pv.Plotter(off_screen=True)  # Close screen display
    # plotter.add_mesh(cavity.mesh_, show_edges=True)
    # plotter.set_background("white")  # Set background colour.
    # plotter.screenshot(path.join(WORK_DIR, "cavity_8_mesh_plot.png"))   # Save figure.
    # plotter.close()

def test_read_cavity_from_xml():
    """Test the `interface.mole2._read_cavity_from_xml` function."""
    xml_filepath = path.join(DATA_DIR, "cavity_calc", "cavities.xml")
    cavity_id = 3  # Use the 3rd cavity for test.
    
    boundary_residues, inner_residues = interface.mole2._read_cavity_from_xml(xml_filepath, cavity_id)
    
    # Verify value type.
    assert isinstance(boundary_residues, str)
    assert isinstance(inner_residues, str)
    
    # Verify content.
    expected_boundary = "ILE 55 A,ARG 58 A,GLY 63 A,PRO 64 A,ARG 126 A,ASP 127 A,ILE 130 A"
    expected_inner = "LEU 50 A,LEU 54 A,ILE 61 A,GLY 62 A,LEU 92 A,SER 94 A,PHE 95 A,TYR 97 A,VAL 125 A,GLY 131 A,GLN 132 A,VAL 160 A,TYR 162 A,ILE 195 A"
    
    assert boundary_residues == expected_boundary
    assert inner_residues == expected_inner
    
    # Test nonexistent cavity_id.
    try:
        interface.mole2._read_cavity_from_xml(xml_filepath, 999)
        assert False, "Should raise ValueError for non-existent cavity_id"
    except ValueError:
        pass

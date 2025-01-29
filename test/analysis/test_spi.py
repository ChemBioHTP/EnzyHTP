"""Testing enzy_htp.analysis.spi.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu> 
Date: 2024-02-24
"""

import glob
import pytest
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_region import create_region_from_selection_pattern
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.analysis import bond_dipole
from enzy_htp import interface
from enzy_htp import PDBParser

from enzy_htp.analysis import spi_metric

from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.core.general import EnablePropagate, get_itself
from enzy_htp.structure.structure_operation import remove_solvent

from enzy_htp._interface.amber_interface import AmberNCParser, AmberRSTParser, AmberMDCRDParser
from enzy_htp.structure.structure_io import PrmtopParser, PDBParser


DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_spi_ligand_sele():
    """To validate that SPI calculation support ligand_sele pattern, 
    which enables SPI calculation can be performed on a substrate or ligand with multiple residues."""
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    prmtop = f"{DATA_DIR}/test_spi.prmtop"
    data = f"{DATA_DIR}/test_spi.mdcrd"
    target_spi:float = 1.555860632049337
    stru = sp.get_structure(pdb_file)

    for tt in stru.residues:
        if tt.name == 'H5J':    
            break

    mdcrd_parser=AmberMDCRDParser(prmtop).get_coordinates

    stru_esm = StructureEnsemble(
        topology=stru,
        top_parser=get_itself,
        coord_parser=mdcrd_parser,
        coordinate_list=data
    )
    
    for single_stru in stru_esm.structures():
        remove_solvent(single_stru)
    
    pocket_sele = "resi 9+11+48+50+101+128+201+202+222"
    ligand_sele = "resn H5J"

    spis = spi_metric(stru_esm, ligand_sele, pocket_sele)

    assert abs(np.mean(np.array(spis)) - target_spi) <= 0.01

def test_spi_consistent_with_old_enzyhtp():
    """To validate that the SPI calculated here is accurate, we compare a value calculated with spi_metric() to a value from the original/old EnzyHTP"""
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    prmtop = f"{DATA_DIR}/test_spi.prmtop"
    data = f"{DATA_DIR}/test_spi.mdcrd"
    target_spi:float = 1.555860632049337

    stru = sp.get_structure(pdb_file)

    for tt in stru.residues:
        if tt.name == 'H5J':    
            break

    mdcrd_parser=AmberMDCRDParser(prmtop).get_coordinates

    se = StructureEnsemble(
        topology=stru,
        top_parser=get_itself,
        coord_parser=mdcrd_parser,
        coordinate_list=data
    )
    
    for ss in se.structures():
        remove_solvent( ss )
    
    aas = set()
    for ii in map( lambda ll: int(ll), "9,11,48,50,101,128,201,202,222".split(',')):
        aas.add(ii)

    active_site = list()
    
    for rr in stru.residues:
        if rr.idx in aas:   
            active_site.append(rr)

    spis = spi_metric( se , tt, active_site)

    assert abs(np.mean(np.array(spis)) - target_spi) <= 0.01


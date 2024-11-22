"""Testing enzy_htp.analysis.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Created: 2024-08-21
"""

# Here put the import lib.
import glob
from typing import List
import pytest
import os
import numpy as np

import enzy_htp.core.file_system as fs
from enzy_htp import PDBParser, interface, _LOGGER

from enzy_htp.analysis import rmsd_calculation, rmsd_of_region, rmsd_of_structure

from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.core.general import EnablePropagate, get_itself

from enzy_htp._interface.amber_interface import AmberNCParser, AmberRSTParser, AmberMDCRDParser
from enzy_htp.structure.structure_io import PrmtopParser, PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

prmtop_path=f"{DATA_DIR}/test_spi.prmtop",
traj_path=f"{DATA_DIR}/test_spi.mdcrd",
ref_pdb=f"{DATA_DIR}/test_spi.pdb"
mdcrd_data = f"{DATA_DIR}/test_spi.mdcrd"
mdcrd_parser=AmberMDCRDParser(prmtop_path).get_coordinates

def test_rmsd_calculation():
    """Test the `rmsd_calculation` API."""
    # structure_ensemble = interface.amber.load_traj(
    #     prmtop_path=prmtop_path,
    #     traj_path=traj_path,
    #     ref_pdb=ref_pdb,
    # )
    structure_ensemble = StructureEnsemble(
        topology=ref_pdb,
        top_parser=sp.get_structure,
        coord_parser=mdcrd_parser,
        coordinate_list=mdcrd_data,
    )
    pattern = "(br. resi 254 around 4 & polymer.protein) and (not elem H)"
    # this equals to 
    # :9,50,101,103,128,130,144,201-202&!@H=
    # :9,11,48,50,101,128,169,201,202,222,224&!@H= (Not match)
    
    answer = 0.5944 # Calculate the rmsd values manually using cpptraj.
    result = rmsd_calculation(
        stru_esm = structure_ensemble,
        region_pattern = pattern,
    )
    # for r, a in zip(result, answer):
    #     assert np.isclose(r, a, atol=0.001)
    assert np.isclose(result, answer, atol=0.001)

def test_rmsd_of_region():
    """Test the `rmsd_of_region` API."""
    structure_ensemble = StructureEnsemble(
        topology=ref_pdb,
        top_parser=sp.get_structure,
        coord_parser=mdcrd_parser,
        coordinate_list=mdcrd_data,
    )
    structure_0 = structure_ensemble.structure_0
    mask_region = list(filter(lambda resi: resi.name.lower()=="ala" , structure_0.residues))

    result = rmsd_of_region(stru_esm=structure_ensemble, mask_region=mask_region, ca_only=False)
    _LOGGER.info(f"The RMSD value is {result}.")

    answer = 0.6252
    assert np.isclose(result, answer, atol=0.001)

def test_rmsd_of_structure():
    """Test the `rmsd_of_structure` API."""
    structure_ensemble = StructureEnsemble(
        topology=ref_pdb,
        top_parser=sp.get_structure,
        coord_parser=mdcrd_parser,
        coordinate_list=mdcrd_data,
    )
    
    result = rmsd_of_structure(stru_esm=structure_ensemble, include_ligand=False, ca_only=False)
    _LOGGER.info(f"The RMSD value is {result}.")

    answer = 1.1462
    assert np.isclose(result, answer, atol=0.001)

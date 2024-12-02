"""Testing enzy_htp.analysis.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Created: 2024-08-21
"""

# Here put the import lib.
from os import path
import glob
import pytest
import numpy as np

import enzy_htp.core.file_system as fs
from enzy_htp import PDBParser, interface, config as eh_config, _LOGGER

from enzy_htp.analysis import rmsd, rmsd_of_structure

from enzy_htp.structure.structure_ensemble import StructureEnsemble

from enzy_htp._interface.amber_interface import AmberMDCRDParser
from enzy_htp.structure.structure_io import PDBParser

DATA_DIR = f"{path.dirname(path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{path.dirname(path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{path.dirname(path.abspath(__file__))}/work_dir/"
sp = PDBParser()

scratch_dir = eh_config["system.SCRATCH_DIR"]
prmtop_path = path.join(DATA_DIR, "test_spi.prmtop")
traj_path = path.join(DATA_DIR, "test_spi.mdcrd")
ref_pdb = path.join(DATA_DIR, "test_spi_chainid.pdb")
mdcrd_data = path.join(DATA_DIR, "test_spi.mdcrd")
mdcrd_parser = AmberMDCRDParser(prmtop_path).get_coordinates

def test_rmsd():
    """Test the `rmsd` API."""
    structure_ensemble: StructureEnsemble = interface.amber.load_traj(
        prmtop_path=prmtop_path,
        traj_path=traj_path,
        ref_pdb=ref_pdb,
    )
    pattern = "(br. resi 254 around 4 & polymer.protein) and (not elem H)"
    # this equals to 
    # :9,50,101,103,128,130,144,201-202&!@H=
    # :9,11,48,50,101,128,169,201,202,222,224&!@H= (Not match)
    
    answer = [  # Calculate the rmsd values manually using cpptraj.
        0.5866,
        0.5253,
        0.6249,
        0.4517,
        0.8457,
        0.5998,
        0.6588,
        0.6615,
        0.5346,
        0.6015,
    ]
    result = rmsd(
        stru_esm=structure_ensemble,
        region_pattern=pattern,
    )
    fs.safe_rmdir(scratch_dir)
    for r, a in zip(result[:len(answer)], answer):
        assert np.isclose(r, a, atol=0.001)

def test_rmsd_of_structure():
    """Test the `rmsd_of_structure` API."""
    structure_ensemble: StructureEnsemble = interface.amber.load_traj(
        prmtop_path=prmtop_path,
        traj_path=traj_path,
        ref_pdb=ref_pdb,
    )
    
    result = rmsd_of_structure(stru_esm=structure_ensemble, include_ligand=False, ca_only=False)
    # _LOGGER.info(f"The RMSD value is {result}.")

    answer = [
        1.1320,
        1.2778,
        1.2113,
        1.1376,
        1.3825,
        1.2082,
        1.1114,
        1.2439,
        1.2094,
        1.2377,
    ]
    fs.safe_rmdir(scratch_dir)
    for r, a in zip(result[:len(answer)], answer):
        assert np.isclose(r, a, atol=0.001)

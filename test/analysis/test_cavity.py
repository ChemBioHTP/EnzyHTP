"""Testing enzy_htp.analysis.cavity

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Created: 2025-06-02
"""

# Here put the import lib.
from os import path
from statistics import mean
import glob
import pytest
import numpy as np

from enzy_htp import interface
from enzy_htp.structure import PDBParser
from enzy_htp.analysis.cavity import identify_stru_cavities, ensemble_cavity_volumes
import enzy_htp.core.file_system as fs

DATA_DIR = f"{path.dirname(path.abspath(__file__))}/data/"
WORK_DIR = f"{path.dirname(path.abspath(__file__))}/work_dir/"
SCRATCH_DIR = f"{path.dirname(path.abspath(__file__))}/scratch/"

amber_interface = interface.amber
sp = PDBParser()

def test_identify_stru_cavities():
    """Test `identify_stru_cavities` function."""
    fs.safe_mkdir(WORK_DIR)
    
    pdb_filepath = path.join(DATA_DIR, "cavity_calc", "aclHMT-ETI-SAH_no-ETI.pdb")
    stru = sp.get_structure(pdb_filepath)
    cavities = identify_stru_cavities(stru=stru, work_dir=WORK_DIR)
    fs.safe_rmdir(WORK_DIR)

    assert len(cavities) == 10
    cavity = cavities[0]
    assert len(cavity.boundary_residues) == 4
    assert len(cavity.inner_residues) == 31
    assert abs(cavity.volume - 1147) < 1
    assert abs(cavity.software_report_volume - 947) < 1

def test_ensemble_cavity_volumes():
    """Test `ensemble_cavity_volumes` function."""
    prmtop_path = path.join(DATA_DIR, "test_spi.prmtop")
    traj_path = path.join(DATA_DIR, "test_spi.mdcrd")
    ref_pdb = path.join(DATA_DIR, "test_spi_chainid.pdb")

    stru_esm = amber_interface.load_traj(
        prmtop_path=prmtop_path,
        traj_path=traj_path,
        ref_pdb=ref_pdb,
    )

    with pytest.raises(ValueError) as exe:
        volumes = ensemble_cavity_volumes(
            stru_esm=stru_esm,
            contain_ligand="resn H5J",
            frame_0_based=True,
            work_dir=SCRATCH_DIR
        )
        assert exe.value
    fs.safe_rmdir(SCRATCH_DIR)
    pass
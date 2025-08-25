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
from enzy_htp import config as eh_config
from enzy_htp.structure import PDBParser
from enzy_htp.analysis.cavity import identify_stru_cavities, ensemble_cavity_volumes
import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_selection import select_stru

DATA_DIR = f"{path.dirname(path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{path.dirname(path.abspath(__file__))}/../structure/data/"
WORK_DIR = f"{path.dirname(path.abspath(__file__))}/work_dir/"
SCRATCH_DIR = f"{path.dirname(path.abspath(__file__))}/scratch/"

amber_interface = interface.amber
sp = PDBParser()

@pytest.fixture
def sample_structure_ensemble(patch_scratch_dir):
    """Fixture to provide a sample StructureEnsemble for testing."""
    stru_esm = interface.amber.load_traj(
        prmtop_path=f"{STRU_DATA_DIR}/KE_07_R7_2_S_10f.prmtop",
        traj_path=f"{STRU_DATA_DIR}/KE_07_R7_2_S_10f.mdcrd",
        ref_pdb=f"{STRU_DATA_DIR}/KE_07_R7_2_S_geom_1.pdb"
    )
    return stru_esm

@pytest.fixture
def patch_scratch_dir(monkeypatch, tmp_path):
    """Fixture to patch the SCRATCH_DIR to a temporary directory for the duration of a test."""
    temp_scratch = tmp_path / "scratch"
    temp_scratch.mkdir()
    monkeypatch.setattr(eh_config.system, 'SCRATCH_DIR', str(temp_scratch))
    yield str(temp_scratch)
    # No need for explicit cleanup, tmp_path and monkeypatch handle it automatically

def test_identify_stru_cavities():
    """Test `identify_stru_cavities` function."""    
    pdb_filepath = path.join(DATA_DIR, "cavity_calc", "aclHMT-ETI-SAH_no-ETI.pdb")
    stru = sp.get_structure(pdb_filepath)
    cavities = identify_stru_cavities(stru=stru, work_dir=WORK_DIR)

    cavity = cavities[0]
    assert len(cavity.boundary_residues) == 4
    assert len(cavity.inner_residues) == 31
    assert abs(cavity.volume - 1147) < 1
    assert abs(cavity.software_report_volume - 947) < 1

def test_ensemble_cavity_volumes(sample_structure_ensemble, patch_scratch_dir):
    """Test `ensemble_cavity_volumes` function."""
    work_dir = f"{WORK_DIR}/cavity"
    stru_esm = sample_structure_ensemble

    volumes = ensemble_cavity_volumes(
        stru_esm=stru_esm,
        contain_ligand="resn H5J",
        frame_0_based=True,
        work_dir=work_dir
    )

    assert len(volumes) == 11
    assert all(isinstance(v, float) for v in volumes)
    assert mean(volumes) > 0


def test_ensemble_cavity_volumes_with_composing_residues(sample_structure_ensemble, patch_scratch_dir):
    """Test `ensemble_cavity_volumes` function with `composing_residues`."""
    work_dir = f"{WORK_DIR}/cavity_composing_residues"

    stru_esm = sample_structure_ensemble

    # Select residues around the ligand in the first frame to define the cavity
    structure_0 = stru_esm.structure_0
    ligand_resn = "H5J"
    composing_residues = select_stru(structure_0, f"byres (resn {ligand_resn} around 5)").involved_residues

    volumes = ensemble_cavity_volumes(
        stru_esm=stru_esm,
        composing_residues=composing_residues,
        frame_0_based=True,
        work_dir=work_dir
    )

    assert len(volumes) == 11
    assert all(isinstance(v, float) for v in volumes)
    assert mean(volumes) > 0


def test_ensemble_cavity_volumes_with_target_cavity(sample_structure_ensemble, patch_scratch_dir):
    """Test `ensemble_cavity_volumes` function with `target_cavity`."""
    work_dir = f"{WORK_DIR}/cavity_target_cavity"

    stru_esm = sample_structure_ensemble

    # Identify cavities in the first frame and select one as the target
    structure_0 = stru_esm.structure_0
    cavities = identify_stru_cavities(stru=structure_0, work_dir=f"{work_dir}/frame_0_cavities")
    target_cavity = cavities[0]

    volumes = ensemble_cavity_volumes(
        stru_esm=stru_esm,
        target_cavity=target_cavity,
        frame_0_based=False,  # Target cavity is already identified
        work_dir=work_dir
    )

    assert len(volumes) == 11
    assert all(isinstance(v, float) for v in volumes)
    assert mean(volumes) > 0


def test_ensemble_cavity_volumes_with_contain_ligand(sample_structure_ensemble, patch_scratch_dir):
    """Test `ensemble_cavity_volumes` function with `contain_ligand`."""
    work_dir = f"{WORK_DIR}/cavity_contain_ligand"

    stru_esm = sample_structure_ensemble

    volumes = ensemble_cavity_volumes(
        stru_esm=stru_esm,
        contain_ligand="resn H5J",
        frame_0_based=True,
        work_dir=work_dir
    )

    assert len(volumes) == 11
    assert all(isinstance(v, float) for v in volumes)
    assert mean(volumes) > 0

"""Testing enzy_htp.analysis.residue_pka.py

Author: QZ Shao <shaoqz@icloud.com>
Date: 2025-08-07
"""
import pytest
import os
import tempfile
import shutil

from enzy_htp import interface
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.analysis import residue_pka
from enzy_htp.core import file_system as fs

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
INT_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../_interface/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"

@pytest.fixture
def sample_structure():
    """Fixture to provide a sample protein structure for testing."""
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    assert os.path.exists(pdb_file), f"Test PDB file not found: {pdb_file}"
    return PDBParser().get_structure(pdb_file)

@pytest.fixture
def sample_structure_ensemble():
    """Fixture to provide a sample StructureEnsemble for testing."""
    # Load trajectory; missing files will raise and cause test failure
    stru_esm = interface.amber.load_traj(
        prmtop_path=f"{INT_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{INT_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{INT_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
    return stru_esm

@pytest.fixture
def temp_work_dir():
    """Fixture to provide a temporary working directory for tests."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir, ignore_errors=True)

def test_residue_pka_structure(sample_structure, temp_work_dir):
    """Test residue pKa calculation for a single Structure."""
    result = residue_pka(
        stru=sample_structure,
        method="propka",
        work_dir=temp_work_dir
    )
    
    # Basic checks
    assert isinstance(result, dict), "Result should be a dictionary"
    assert len(result) > 0, "Result should contain pKa values"

    # Check that all keys are (chain_id, res_num) tuples and values are numbers
    for key, pka_val in result.items():
        assert isinstance(key, tuple), f"Key {key} should be a tuple"
        assert len(key) == 2, f"Key {key} should be a 2-tuple"
        chain_id, res_num = key
        assert isinstance(chain_id, str), f"Chain ID {chain_id} should be a string"
        assert isinstance(res_num, int), f"Residue number {res_num} should be an integer"
        assert isinstance(pka_val, (int, float)), f"pKa value {pka_val} should be numeric"
        assert 0 < pka_val < 20, f"pKa value {pka_val} should be reasonable (0-20)"

def test_residue_pka_structure_with_target_residues(sample_structure, temp_work_dir):
    """Test residue pKa calculation for specific target residues."""
    # First get all pKa values to see what's available
    all_pka = residue_pka(
        stru=sample_structure,
        method="propka",
        work_dir=temp_work_dir
    )
    
    if len(all_pka) == 0:
        pytest.skip("No ionizable residues found in test structure")
    
    # Test with specific residues using chain.residue format
    target_keys = list(all_pka.keys())[:2]  # Take first two residues
    target_residues = [f"{chain_id}.{res_num}" for chain_id, res_num in target_keys]
    result = residue_pka(
        stru=sample_structure,
        target_residues=target_residues,
        method="propka",
        work_dir=temp_work_dir
    )
    
    assert len(result) <= len(target_keys), "Result should not exceed requested residues"
    for key in result.keys():
        assert key in target_keys, f"Residue {key} not in target list"

def test_residue_pka_structure_ensemble(sample_structure_ensemble, temp_work_dir):
    """Test residue pKa calculation for a StructureEnsemble."""
    result = residue_pka(
        stru=sample_structure_ensemble,
        method="propka",
        work_dir=temp_work_dir
    )
    
    # Basic checks
    assert isinstance(result, list), "Result should be a list for StructureEnsemble"
    assert len(result) > 0, "Result should contain pKa dictionaries for each structure"
    
    # Check each structure's results
    for i, struct_result in enumerate(result):
        assert isinstance(struct_result, dict), f"Structure {i} result should be a dictionary"
        for key, pka_val in struct_result.items():
            assert isinstance(key, tuple), f"Key {key} should be a tuple"
            assert len(key) == 2, f"Key {key} should be a 2-tuple"
            chain_id, res_num = key
            assert isinstance(chain_id, str), f"Chain ID {chain_id} should be a string"
            assert isinstance(res_num, int), f"Residue number {res_num} should be an integer"
            assert isinstance(pka_val, (int, float)), f"pKa value {pka_val} should be numeric"
    
    fs.clean_temp_file_n_dir([
        sample_structure_ensemble.topology_source_file,
    ])

def test_residue_pka_invalid_method(sample_structure, temp_work_dir):
    """Test that invalid methods raise appropriate errors."""
    with pytest.raises(ValueError):
        residue_pka(
            stru=sample_structure,
            method="invalid_method",
            work_dir=temp_work_dir
        )

def test_residue_pka_invalid_structure_type(temp_work_dir):
    """Test that invalid structure types raise appropriate errors."""
    with pytest.raises(TypeError):
        residue_pka(
            stru="not_a_structure",
            method="propka",
            work_dir=temp_work_dir
        )

def test_propka_interface_directly(sample_structure, temp_work_dir):
    """Test the PROPKA interface directly to ensure it works."""
    # Remove solvents first to avoid PDB formatting issues
    from enzy_htp.preparation.clean import remove_solvent
    clean_structure = remove_solvent(sample_structure, in_place=False)
    
    result = interface.propka.get_residue_pka_from_stru(
        stru=clean_structure,
        work_dir=temp_work_dir
    )
    
    assert isinstance(result, dict), "PROPKA interface should return a dictionary"
    # Note: We don't assert on specific values as they depend on the test structure
    # The user mentioned they will manually verify correct pKa results later

def test_propka_methods_registry():
    """Test that the PKA_METHODS registry is properly set up."""
    from enzy_htp.analysis.residue_pka import PKA_METHODS
    
    assert "propka" in PKA_METHODS, "PROPKA method should be registered"
    assert callable(PKA_METHODS["propka"]), "PROPKA method should be callable"
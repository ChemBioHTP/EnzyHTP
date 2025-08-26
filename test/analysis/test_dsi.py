"""Test module for enzy_htp.analysis.dsi

Author: QZ Shao <shaoqz@icloud.com>
Date: 2025-08-15
"""
import os
import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from enzy_htp.analysis.dsi import dsi, _expand_residue_ranges
from enzy_htp.structure import StructureEnsemble, PDBParser
from enzy_htp import interface
from enzy_htp import config as eh_config

# Test data directory
DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"

@pytest.fixture
def patch_scratch_dir(monkeypatch, tmp_path):
    """Fixture to patch the SCRATCH_DIR to a temporary directory for the duration of a test."""
    temp_scratch = tmp_path / "scratch"
    temp_scratch.mkdir()
    monkeypatch.setattr(eh_config.system, 'SCRATCH_DIR', str(temp_scratch))
    yield str(temp_scratch)
    # No need for explicit cleanup, tmp_path and monkeypatch handle it automatically

def test_expand_residue_ranges():
    """Test expansion of residue ranges."""
    # Test empty list
    assert _expand_residue_ranges([]) == []
    
    # Test single residue
    result = _expand_residue_ranges([("A", 10)])
    assert result == [("A", 10)]
    
    # Test two residues from same chain - should expand to range
    result = _expand_residue_ranges([("A", 10), ("A", 15)])
    expected = [("A", 10), ("A", 11), ("A", 12), ("A", 13), ("A", 14), ("A", 15)]
    assert result == expected
    
    # Test two residues from different chains - should remain as individual residues
    result = _expand_residue_ranges([("A", 10), ("B", 15)])
    expected = [("A", 10), ("B", 15)]
    assert result == expected
    
    # Test multiple residues - should use individual residues (no expansion)
    result = _expand_residue_ranges([("A", 10), ("A", 12), ("A", 15)])
    expected = [("A", 10), ("A", 12), ("A", 15)]
    assert result == expected


def test_dsi_input_validation():
    """Test that dsi function validates inputs correctly."""
    # Test with non-StructureEnsemble input
    with pytest.raises(TypeError, match="ensemble must be a StructureEnsemble object"):
        dsi("not_an_ensemble", [("A", 1)], [("A", 2)])
    
    # Create mock ensemble for other tests
    mock_ensemble = MagicMock(spec=StructureEnsemble)
    
    # Test with unsupported engine
    with pytest.raises(ValueError, match="Unsupported engine"):
        dsi(mock_ensemble, [("A", 1)], [("A", 2)], engine="unsupported")
    
    # Test with empty domain residue lists
    with pytest.raises(ValueError, match="Domain residue lists cannot be empty"):
        dsi(mock_ensemble, [], [("A", 2)])
    
    with pytest.raises(ValueError, match="Domain residue lists cannot be empty"):
        dsi(mock_ensemble, [("A", 1)], [])


def test_dsi_function_basic():
    """Test basic functionality of dsi function."""
    with patch('enzy_htp.analysis.dsi.DSI_METHODS') as mock_dsi_methods:
        # Mock the DSI methods dictionary
        expected_result = np.array([1.5, 2.0, 1.8])
        mock_calculate_dsi = MagicMock(return_value=expected_result)
        mock_dsi_methods.__getitem__.return_value = mock_calculate_dsi
        mock_dsi_methods.__contains__.return_value = True
        
        # Create mock ensemble
        mock_ensemble = MagicMock(spec=StructureEnsemble)
        
        # Test with range expansion
        domain1_residues = [("A", 10), ("A", 15)]  # Should expand to range
        domain2_residues = [("B", 20), ("B", 25)]  # Should expand to range

        result = dsi(mock_ensemble, domain1_residues, domain2_residues)
        
        # Verify the mock was called with expanded residues
        mock_calculate_dsi.assert_called_once()
        args = mock_calculate_dsi.call_args[0]
        
        assert args[0] == mock_ensemble
        # Check that residues were expanded
        assert len(args[1]) == 6  # 10-15 inclusive
        assert len(args[2]) == 6  # 20-25 inclusive
        
        # Check return value
        np.testing.assert_array_equal(result, expected_result)


def test_dsi_function_individual_residues():
    """Test dsi function with individual residues (no range expansion)."""
    with patch('enzy_htp.analysis.dsi.DSI_METHODS') as mock_dsi_methods:
        # Mock the DSI methods dictionary
        expected_result = np.array([3.2, 2.8, 3.1])
        mock_calculate_dsi = MagicMock(return_value=expected_result)
        mock_dsi_methods.__getitem__.return_value = mock_calculate_dsi
        mock_dsi_methods.__contains__.return_value = True
        
        # Create mock ensemble
        mock_ensemble = MagicMock(spec=StructureEnsemble)
        
        # Test with individual residues (more than 2)
        domain1_residues = [("A", 10), ("A", 15), ("A", 20)]
        domain2_residues = [("B", 30), ("B", 35), ("B", 40)]

        result = dsi(mock_ensemble, domain1_residues, domain2_residues)
        
        # Verify the mock was called with the same residues (no expansion)
        mock_calculate_dsi.assert_called_once()
        args = mock_calculate_dsi.call_args[0]
        
        assert args[0] == mock_ensemble
        assert args[1] == domain1_residues  # Should be unchanged
        assert args[2] == domain2_residues  # Should be unchanged
        
        # Check return value
        np.testing.assert_array_equal(result, expected_result)


def test_dsi_function_default_engine():
    """Test that dsi function uses cpptraj as default engine."""
    with patch('enzy_htp.analysis.dsi.DSI_METHODS') as mock_dsi_methods:
        # Mock the DSI methods dictionary
        expected_result = np.array([1.0])
        mock_calculate_dsi = MagicMock(return_value=expected_result)
        mock_dsi_methods.__getitem__.return_value = mock_calculate_dsi
        mock_dsi_methods.__contains__.return_value = True
        
        # Create mock ensemble
        mock_ensemble = MagicMock(spec=StructureEnsemble)
        
        # Test with default engine (should be cpptraj)
        domain1_residues = [("A", 1)]
        domain2_residues = [("A", 2)]

        result = dsi(mock_ensemble, domain1_residues, domain2_residues)
        
        # Verify the DSI_METHODS dictionary was accessed with "cpptraj"
        mock_dsi_methods.__getitem__.assert_called_once_with("cpptraj")
        
        # Check return value
        np.testing.assert_array_equal(result, expected_result)


def test_dsi_methods_dictionary():
    """Test that DSI_METHODS contains expected engines."""
    from enzy_htp.analysis.dsi import DSI_METHODS
    from enzy_htp import interface
    
    assert "cpptraj" in DSI_METHODS
    assert DSI_METHODS["cpptraj"] == interface.amber.calculate_dsi_metrics


def test_dsi_real_data(patch_scratch_dir):
    """Test DSI calculation with real trajectory data."""
    # Use the same test data as other analysis tests
    prmtop_path = os.path.join(DATA_DIR, "test_spi.prmtop")
    traj_path = os.path.join(DATA_DIR, "test_spi.mdcrd")
    ref_pdb = os.path.join(DATA_DIR, "test_spi_chainid.pdb")
    
    # Load trajectory ensemble
    structure_ensemble = interface.amber.load_traj(
        prmtop_path=prmtop_path,
        traj_path=traj_path,
        ref_pdb=ref_pdb,
    )
    
    # Define two domains for DSI calculation
    # Domain 1: residues 1-10 (N-terminal region)
    # Domain 2: residues 100-110 (middle region)
    domain1_residues = [("A", 1), ("A", 10)]
    domain2_residues = [("A", 100), ("A", 110)]
    
    # Calculate DSI using the analysis API
    result = dsi(structure_ensemble, domain1_residues, domain2_residues)

    # Validate results
    assert isinstance(result, np.ndarray)
    assert len(result) > 0  # Should have at least one frame
    assert all(isinstance(x, (int, float)) for x in result)  # All values should be numeric
    
    # DSI values should be reasonable (not NaN or infinite)
    assert not np.any(np.isnan(result))
    assert not np.any(np.isinf(result))
    
    # DSI values should generally be positive (distance minus radii of gyration)
    # but can be negative if domains overlap significantly
    assert all(x > -50 for x in result)  # Reasonable lower bound
    assert all(x < 200 for x in result)   # Reasonable upper bound

def test_dsi_real_data_single_stru(patch_scratch_dir):
    """Test DSI calculation with real trajectory data."""
    sp = PDBParser()
    # Use the same test data as other analysis tests
    ref_pdb = os.path.join(DATA_DIR, "test_spi_chainid.pdb")
    
    # Load trajectory ensemble
    structure_ensemble = StructureEnsemble(
        topology = ref_pdb,
        top_parser = sp.get_structure,
        coordinate_list = ref_pdb,
        coord_parser = sp.get_structure,
    )
    
    # Define two domains for DSI calculation
    # Domain 1: residues 1-10 (N-terminal region)
    # Domain 2: residues 100-110 (middle region)
    domain1_residues = [("A", 1), ("A", 10)]
    domain2_residues = [("A", 100), ("A", 110)]
    
    # Calculate DSI using the analysis API
    result = dsi(structure_ensemble, domain1_residues, domain2_residues)

    # Validate results
    assert isinstance(result, np.ndarray)
    assert len(result) == 1  # Should have exactly one frame
    assert all(isinstance(x, (int, float)) for x in result)  # All values should be numeric
    
    # DSI values should be reasonable (not NaN or infinite)
    assert not np.any(np.isnan(result))
    assert not np.any(np.isinf(result))
    
    # DSI values should generally be positive (distance minus radii of gyration)
    # but can be negative if domains overlap significantly
    assert all(x > -50 for x in result)  # Reasonable lower bound
    assert all(x < 200 for x in result)   # Reasonable upper bound


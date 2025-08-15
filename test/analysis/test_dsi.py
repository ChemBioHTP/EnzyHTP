"""Test module for enzy_htp.analysis.dsi

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-08-15
"""
import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from enzy_htp.analysis.dsi import dsi, _expand_residue_ranges
from enzy_htp.structure import StructureEnsemble


def test_expand_residue_ranges():
    """Test expansion of residue ranges."""
    # Test empty list
    assert _expand_residue_ranges([]) == []
    
    # Test single residue
    result = _expand_residue_ranges([("A", 10)])
    assert result == [("A", 10)]
    
    # Test two residues - should expand to range
    result = _expand_residue_ranges([("A", 10), ("A", 15)])
    expected = [("A", 10), ("A", 11), ("A", 12), ("A", 13), ("A", 14), ("A", 15)]
    assert result == expected
    
    # Test multiple residues - should use individual residues
    result = _expand_residue_ranges([("A", 10), ("A", 12), ("A", 15)])
    expected = [("A", 10), ("A", 12), ("A", 15)]
    assert result == expected
    
    # Test multiple chains with ranges
    result = _expand_residue_ranges([("A", 10), ("A", 12), ("B", 20), ("B", 22)])
    expected = [("A", 10), ("A", 12), ("B", 20), ("B", 22)]
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
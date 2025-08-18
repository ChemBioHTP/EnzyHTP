"""Test module for enzy_htp.structure_prediction.prediction

Author: Gemini
Date: 2025-08-16
"""
import pytest
from unittest.mock import MagicMock, patch, Mock
from typing import Dict

from enzy_htp.structure_prediction.prediction import predict_structure, PREDICTION_ENGINES
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJobConfig


def test_predict_structure_unsupported_engine():
    """Test predict_structure with an unsupported engine."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    # Act & Assert
    with pytest.raises(ValueError, match="Unsupported engine: unsupported_engine"):
        predict_structure(sequence, engine="unsupported_engine")


@patch('enzy_htp._interface.alphafold_interface.subprocess.run')
@patch('enzy_htp.structure_prediction.prediction.af2_predict')
def test_predict_structure_single_sequence(mock_af2_predict, mock_subprocess):
    """Test predict_structure with a single sequence."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    mock_structure = Mock(spec=Structure)
    mock_af2_predict.return_value = {"seq_0": mock_structure}
    
    # Act
    result = predict_structure(sequence, engine="alphafold2")
    
    # Assert
    assert isinstance(result, dict)
    assert "seq_0" in result
    assert result["seq_0"] == mock_structure
    mock_af2_predict.assert_called_once_with([sequence], None)


@patch('enzy_htp.structure_prediction.prediction.af2_predict')
def test_predict_structure_multiple_sequences(mock_af2_predict):
    """Test predict_structure with multiple sequences."""
    # Arrange
    sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY"]
    mock_structures = {
        "seq_0": Mock(spec=Structure),
        "seq_1": Mock(spec=Structure)
    }
    mock_af2_predict.return_value = mock_structures
    
    # Act
    result = predict_structure(sequences, engine="alphafold2")
    
    # Assert
    assert isinstance(result, dict)
    assert len(result) == 2
    assert "seq_0" in result
    assert "seq_1" in result
    mock_af2_predict.assert_called_once_with(sequences, None)


@patch('enzy_htp.structure_prediction.prediction.af2_predict')
def test_predict_structure_with_cluster_job_config(mock_af2_predict):
    """Test predict_structure with cluster job configuration."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    cluster_config = ClusterJobConfig(res_keywords={"partition": "gpu", "walltime": "02:00:00"})
    mock_structure = Mock(spec=Structure)
    mock_af2_predict.return_value = {"seq_0": mock_structure}
    
    # Act
    result = predict_structure(sequence, engine="alphafold2", cluster_job_config=cluster_config)
    
    # Assert
    assert isinstance(result, dict)
    mock_af2_predict.assert_called_once_with([sequence], cluster_config)


@patch('enzy_htp.structure_prediction.prediction.af2_predict')
def test_predict_structure_with_kwargs(mock_af2_predict):
    """Test predict_structure with additional keyword arguments."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    mock_structure = Mock(spec=Structure)
    mock_af2_predict.return_value = {"seq_0": mock_structure}
    
    # Act
    result = predict_structure(
        sequence, 
        engine="alphafold2", 
        num_models=3, 
        num_recycles=5
    )
    
    # Assert
    assert isinstance(result, dict)
    mock_af2_predict.assert_called_once_with(
        [sequence], None, num_models=3, num_recycles=5
    )


def test_prediction_engines_registry():
    """Test that the PREDICTION_ENGINES registry is properly configured."""
    # Assert
    assert "alphafold2" in PREDICTION_ENGINES
    assert callable(PREDICTION_ENGINES["alphafold2"])

"""Test module for enzy_htp.structure_prediction.prediction

Author: QZ Shao <shaoqz@icloud.com>
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


@patch('enzy_htp.structure_prediction.prediction.PREDICTION_ENGINES')
def test_predict_structure_single_sequence(mock_engines):
    """Test predict_structure with a single sequence."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    mock_structure = Mock(spec=Structure)
    mock_af2_predict = Mock(return_value={sequence: mock_structure})
    mock_engines.__getitem__.return_value = mock_af2_predict
    mock_engines.__contains__.return_value = True
    
    # Act
    result = predict_structure(sequence, engine="alphafold2")
    
    # Assert
    assert isinstance(result, dict)
    assert sequence in result
    assert result[sequence] == mock_structure
    mock_af2_predict.assert_called_once_with([sequence], cluster_job_config=None)


@patch('enzy_htp.structure_prediction.prediction.PREDICTION_ENGINES')
def test_predict_structure_multiple_sequences(mock_engines):
    """Test predict_structure with multiple sequences."""
    # Arrange
    sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY"]
    mock_structures = {
        sequences[0]: Mock(spec=Structure),
        sequences[1]: Mock(spec=Structure)
    }
    mock_af2_predict = Mock(return_value=mock_structures)
    mock_engines.__getitem__.return_value = mock_af2_predict
    mock_engines.__contains__.return_value = True
    
    # Act
    result = predict_structure(sequences, engine="alphafold2")
    
    # Assert
    assert isinstance(result, dict)
    assert len(result) == 2
    assert sequences[0] in result
    assert sequences[1] in result
    mock_af2_predict.assert_called_once_with(sequences, cluster_job_config=None)


@patch('enzy_htp.structure_prediction.prediction.PREDICTION_ENGINES')
def test_predict_structure_with_cluster_job_config(mock_engines):
    """Test predict_structure with cluster job configuration."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    cluster_config = ClusterJobConfig(res_keywords={"partition": "gpu", "walltime": "02:00:00"})
    mock_structure = Mock(spec=Structure)
    mock_af2_predict = Mock(return_value={sequence: mock_structure})
    mock_engines.__getitem__.return_value = mock_af2_predict
    mock_engines.__contains__.return_value = True
    
    # Act
    result = predict_structure(sequence, engine="alphafold2", cluster_job_config=cluster_config)
    
    # Assert
    assert isinstance(result, dict)
    mock_af2_predict.assert_called_once_with([sequence], cluster_job_config=cluster_config)


@patch('enzy_htp.structure_prediction.prediction.PREDICTION_ENGINES')
def test_predict_structure_with_kwargs(mock_engines):
    """Test predict_structure with additional keyword arguments."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    mock_structure = Mock(spec=Structure)
    mock_af2_predict = Mock(return_value={sequence: mock_structure})
    mock_engines.__getitem__.return_value = mock_af2_predict
    mock_engines.__contains__.return_value = True
    
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
        [sequence], cluster_job_config=None, num_models=3, num_recycles=5
    )


@patch('enzy_htp.structure_prediction.prediction.PREDICTION_ENGINES')
@patch('enzy_htp.structure_prediction.prediction._parse_sequences_input')
def test_predict_structure_with_fasta_file(mock_parse_input, mock_engines):
    """Test predict_structure with FASTA file input."""
    # Arrange
    from pathlib import Path
    fasta_file = Path("/fake/path/to/sequences.fasta")
    sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY"]
    mock_structures = {
        sequences[0]: Mock(spec=Structure),
        sequences[1]: Mock(spec=Structure)
    }
    
    mock_parse_input.return_value = sequences
    mock_af2_predict = Mock(return_value=mock_structures)
    mock_engines.__getitem__.return_value = mock_af2_predict
    mock_engines.__contains__.return_value = True
    
    # Act
    result = predict_structure(fasta_file, engine="alphafold2")
    
    # Assert
    assert isinstance(result, dict)
    assert len(result) == 2
    mock_parse_input.assert_called_once_with(fasta_file)
    mock_af2_predict.assert_called_once_with(sequences, cluster_job_config=None)


def test_prediction_engines_registry():
    """Test that the PREDICTION_ENGINES registry is properly configured."""
    # Assert
    assert "alphafold2" in PREDICTION_ENGINES
    assert callable(PREDICTION_ENGINES["alphafold2"])

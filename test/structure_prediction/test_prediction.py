"""Test module for enzy_htp.structure_prediction.prediction

Author: Gemini
Date: 2025-08-16
"""
import pytest
from unittest.mock import MagicMock, patch

from enzy_htp.structure_prediction.prediction import predict_structure
from enzy_htp.structure import Structure


@patch('enzy_htp.interface.alphafold.run')
def test_predict_structure_alphafold(mock_alphafold_run):
    """Test predict_structure with the alphafold engine."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    mock_structure = MagicMock(spec=Structure)
    mock_alphafold_run.return_value = mock_structure

    # Act
    result = predict_structure(sequence)

    # Assert
    mock_alphafold_run.assert_called_once_with(sequence)
    assert result is mock_structure


def test_predict_structure_unsupported_engine():
    """Test predict_structure with an unsupported engine."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    # Act & Assert
    with pytest.raises(ValueError, match="Unsupported engine: unsupported_engine"):
        predict_structure(sequence, engine="unsupported_engine")

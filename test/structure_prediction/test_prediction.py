"""Test module for enzy_htp.structure_prediction.prediction

Author: Gemini
Date: 2025-08-16
"""
import pytest
from unittest.mock import MagicMock, patch

from enzy_htp.structure_prediction.prediction import predict_structure
from enzy_htp.structure import Structure

def test_predict_structure_unsupported_engine():
    """Test predict_structure with an unsupported engine."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    # Act & Assert
    with pytest.raises(ValueError, match="Unsupported engine: unsupported_engine"):
        predict_structure(sequence, engine="unsupported_engine")

# TODO make a test without any mock

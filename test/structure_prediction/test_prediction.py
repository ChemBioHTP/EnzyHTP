"""Test module for enzy_htp.structure_prediction.prediction

Author: QZ Shao <shaoqz@icloud.com>
Date: 2025-08-16
"""
import pytest
from unittest.mock import MagicMock, patch, Mock
from typing import Dict
import os
from pathlib import Path

from enzy_htp import interface
from enzy_htp import config as eh_config
from enzy_htp.structure_prediction.prediction import predict_structure
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJobConfig
import enzy_htp.core.file_system as fs
af_interface = interface.alphafold
af_config = eh_config.alphafold
DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"

@pytest.fixture
def alphafold_config_modifier():
    """Fixture to modify any attributes in alphafold config and restore them after test."""
    original_values = {}
    
    def _modify_config(**kwargs):
        """Modify config attributes and store original values for restoration.
        
        Args:
            **kwargs: Key-value pairs where keys are attribute names and values are new values.
            
        Returns:
            The modified config object.
            
        Examples:
            # Modify single attribute
            config = alphafold_config_modifier(INSTALL_TYPE="alphafold2_native_python")
            
            # Modify multiple attributes
            config = alphafold_config_modifier(
                INSTALL_TYPE="alphafold2_native_python",
                EXECUTABLE_PATH="/path/to/run_alphafold.py",
                DATA_DIR="/path/to/data"
            )
        """
        for attr_name, new_value in kwargs.items():
            if hasattr(af_config, attr_name):
                # Store original value if not already stored
                if attr_name not in original_values:
                    original_values[attr_name] = getattr(af_config, attr_name)
                # Set new value
                setattr(af_config, attr_name, new_value)
            else:
                raise AttributeError(f"AlphafoldConfig has no attribute '{attr_name}'")
        return af_config
    
    yield _modify_config
    
    # Restore all original values after test
    for attr_name, original_value in original_values.items():
        setattr(af_config, attr_name, original_value)

def test_predict_structure_unsupported_engine():
    """Test predict_structure with an unsupported engine."""
    # Arrange
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    # Act & Assert
    with pytest.raises(ValueError, match="Unsupported engine: unsupported_engine"):
        predict_structure(sequence, engine="unsupported_engine")

def test_af2_predict_colabfold():
    """Test predict_structure without mocking use colabfold"""
    test_output_dir = Path(f"{WORK_DIR}/test_stru_pred_output")
    # Arrange
    sequence = "MSTPSLIPSGVHEVLAKYKDGN"
    # Act
    result = predict_structure(sequence, work_dir=test_output_dir)

    # Assert
    assert isinstance(result, dict)
    # The result keys should be the actual sequences, not seq_0
    assert sequence in result
    
    # Check comprehensive output structure
    seq_result = result[sequence]
    assert isinstance(seq_result, dict)
    
    # Check best model
    assert "best_model" in seq_result
    assert "best_model_plddt" in seq_result
    assert "best_model_index" in seq_result
    assert seq_result["best_model"].num_residues == 22
    assert seq_result["best_model"].num_atoms == 163
    
    # Check individual models exist
    assert "model_1" in seq_result or "model_2" in seq_result or "model_3" in seq_result or "model_4" in seq_result or "model_5" in seq_result
    
    # Check pLDDT scores are valid
    assert isinstance(seq_result["best_model_plddt"], list)
    assert len(seq_result["best_model_plddt"]) == 22  # Same as number of residues
    fs.safe_rmdir(test_output_dir)


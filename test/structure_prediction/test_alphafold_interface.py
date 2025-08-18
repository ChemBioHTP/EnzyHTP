"""Test module for enzy_htp._interface.alphafold_interface

Author: Claude Code
Date: 2025-08-18
"""
import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import MagicMock, patch, Mock, mock_open

from enzy_htp._interface.alphafold_interface import AlphafoldInterface, af2_predict
from enzy_htp._config.alphafold_config import AlphafoldConfig
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig


class TestAlphafoldInterface:
    """Test class for AlphafoldInterface."""

    def setup_method(self):
        """Set up test fixtures."""
        self.mock_parent = Mock()
        self.config = AlphafoldConfig()
        self.interface = AlphafoldInterface(self.mock_parent, self.config)

    def test_init(self):
        """Test AlphafoldInterface initialization."""
        assert self.interface.config_ == self.config
        assert self.interface.parent_ == self.mock_parent

    @patch('enzy_htp._interface.alphafold_interface.tempfile.NamedTemporaryFile')
    @patch('enzy_htp._interface.alphafold_interface.tempfile.mkdtemp')
    @patch('enzy_htp._interface.alphafold_interface.PDBParser')
    @patch('enzy_htp._interface.alphafold_interface.os.unlink')
    def test_predict_single_sequence_local(self, mock_unlink, mock_pdb_parser, 
                                          mock_mkdtemp, mock_temp_file):
        """Test predict method with single sequence running locally."""
        # Arrange
        sequences = ["ACDEFGHIKLMNPQRSTVWY"]
        mock_temp_file.return_value.__enter__.return_value.name = "/tmp/test.fasta"
        mock_mkdtemp.return_value = "/tmp/alphafold_test"
        
        mock_structure = Mock(spec=Structure)
        mock_parser = Mock()
        mock_parser.read_structure.return_value = mock_structure
        mock_pdb_parser.return_value = mock_parser
        
        # Mock the run method
        self.interface.run = Mock(return_value={"seq_0": "/tmp/alphafold_test/seq_0.pdb"})
        
        # Act
        result = self.interface.predict(sequences)
        
        # Assert
        assert isinstance(result, dict)
        assert "seq_0" in result
        assert result["seq_0"] == mock_structure
        mock_unlink.assert_called_once()

    @patch('enzy_htp._interface.alphafold_interface.tempfile.NamedTemporaryFile')
    @patch('enzy_htp._interface.alphafold_interface.tempfile.mkdtemp')
    def test_predict_with_cluster_job(self, mock_mkdtemp, mock_temp_file):
        """Test predict method with cluster job configuration."""
        # Arrange
        sequences = ["ACDEFGHIKLMNPQRSTVWY"]
        cluster_config = ClusterJobConfig(res_keywords={"partition": "gpu"})
        mock_temp_file.return_value.__enter__.return_value.name = "/tmp/test.fasta"
        mock_mkdtemp.return_value = "/tmp/alphafold_test"
        
        # Mock the make_job method and job behavior
        mock_job = Mock(spec=ClusterJob)
        mock_job.get_state.return_value = ("completed", "completed")
        self.interface.make_job = Mock(return_value=mock_job)
        self.interface._find_output_files = Mock(return_value={"seq_0": "/tmp/test.pdb"})
        
        # Mock PDB parser
        with patch('enzy_htp._interface.alphafold_interface.PDBParser') as mock_pdb_parser:
            mock_structure = Mock(spec=Structure)
            mock_parser = Mock()
            mock_parser.read_structure.return_value = mock_structure
            mock_pdb_parser.return_value = mock_parser
            
            # Act
            result = self.interface.predict(sequences, cluster_job_config=cluster_config)
            
            # Assert
            assert isinstance(result, dict)
            mock_job.submit.assert_called_once()
            mock_job.wait_to_finish.assert_called_once()

    def test_build_colabfold_container_command(self):
        """Test building ColabFold container command."""
        # Arrange
        self.config.INSTALL_TYPE = "colabfold_container"
        self.config.CONTAINER_TYPE = "apptainer"
        self.config.CONTAINER_PATH = "/path/to/container.sif"
        self.config.CONTAINER_BIND_PATHS = {"/cache": "/host/cache", "/work": "/host/work"}
        
        # Act
        cmd = self.interface._build_colabfold_container_command(
            "/tmp/test.fasta", Path("/tmp/out"), 3, 2, 1, 200, False, "alphafold2_ptm", None
        )
        
        # Assert
        assert cmd[0] == "apptainer"
        assert "run" in cmd
        assert "--nv" in cmd
        assert "colabfold_batch" in cmd
        assert "/tmp/test.fasta" in cmd
        assert "--num-models" in cmd
        assert "3" in cmd
        assert "--num-recycle" in cmd
        assert "2" in cmd

    def test_build_alphafold_native_command(self):
        """Test building native AlphaFold command."""
        # Arrange
        self.config.INSTALL_TYPE = "alphafold_native"
        self.config.EXECUTABLE_PATH = "/path/to/run_alphafold.py"
        self.config.DATA_DIR = "/path/to/data"
        
        # Act
        cmd = self.interface._build_alphafold_native_command(
            "/tmp/test.fasta", Path("/tmp/out"), "2023-05-01", "monomer_ptm", "full_dbs", None
        )
        
        # Assert
        assert cmd[0] == "python"
        assert "/path/to/run_alphafold.py" in cmd
        assert "--fasta_paths" in cmd
        assert "/tmp/test.fasta" in cmd
        assert "--data_dir" in cmd
        assert "/path/to/data" in cmd
        assert "--max_template_date" in cmd
        assert "2023-05-01" in cmd

    def test_build_colabfold_python_command(self):
        """Test building ColabFold Python command."""
        # Arrange
        self.config.INSTALL_TYPE = "colabfold_python"
        self.config.EXECUTABLE_PATH = "colabfold_batch"
        
        # Act
        cmd = self.interface._build_colabfold_python_command(
            "/tmp/test.fasta", Path("/tmp/out"), 2, 3, 1, 100, True, "alphafold2", None
        )
        
        # Assert
        assert cmd[0] == "colabfold_batch"
        assert "/tmp/test.fasta" in cmd
        assert "--num-models" in cmd
        assert "2" in cmd
        assert "--templates" in cmd
        assert "--amber" in cmd

    @patch('enzy_htp._interface.alphafold_interface.subprocess.run')
    def test_run_colabfold_container_success(self, mock_subprocess):
        """Test successful run with ColabFold container."""
        # Arrange
        self.config.INSTALL_TYPE = "colabfold_container"
        mock_subprocess.return_value.returncode = 0
        mock_subprocess.return_value.stdout = "Success"
        self.interface._find_output_files = Mock(return_value={"seq_0": "/tmp/out.pdb"})
        
        # Act
        result = self.interface.run("/tmp/test.fasta", "/tmp/out")
        
        # Assert
        assert isinstance(result, dict)
        mock_subprocess.assert_called_once()

    @patch('enzy_htp._interface.alphafold_interface.subprocess.run')
    def test_run_failure(self, mock_subprocess):
        """Test run method with subprocess failure."""
        # Arrange
        from subprocess import CalledProcessError
        mock_subprocess.side_effect = CalledProcessError(1, "cmd", stderr="Error message")
        
        # Act & Assert
        with pytest.raises(RuntimeError, match="AlphaFold execution failed"):
            self.interface.run("/tmp/test.fasta", "/tmp/out")

    def test_make_job(self):
        """Test make_job method."""
        # Arrange
        cluster_config = ClusterJobConfig(res_keywords={"partition": "gpu", "walltime": "02:00:00"})
        self.interface._build_colabfold_container_command = Mock(return_value=["cmd", "args"])
        
        # Act
        with patch('enzy_htp._interface.alphafold_interface.ClusterJob') as mock_cluster_job:
            job = self.interface.make_job(
                "/tmp/test.fasta", "/tmp/out", cluster_job_config=cluster_config
            )
            
            # Assert
            mock_cluster_job.assert_called_once()

    def test_find_output_files(self):
        """Test _find_output_files method."""
        # Arrange
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test PDB files
            (temp_path / "seq_0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            (temp_path / "seq_0_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            (temp_path / "seq_1_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            
            # Act
            result = self.interface._find_output_files(temp_path)
            
            # Assert
            assert isinstance(result, dict)
            assert "seq_0" in result
            assert "seq_1" in result
            # Should prefer relaxed over unrelaxed
            assert "relaxed" in result["seq_0"]

    def test_unsupported_install_type(self):
        """Test run method with unsupported install type."""
        # Arrange
        self.config.INSTALL_TYPE = "unsupported"
        
        # Act & Assert
        with pytest.raises(ValueError, match="Unsupported install type"):
            self.interface.run("/tmp/test.fasta", "/tmp/out")


class TestAf2Predict:
    """Test class for af2_predict function."""

    @patch('enzy_htp._interface.alphafold_interface.interface.alphafold.predict')
    def test_af2_predict_single_sequence(self, mock_predict):
        """Test af2_predict with single sequence."""
        # Arrange
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        mock_structure = Mock(spec=Structure)
        mock_predict.return_value = {"seq_0": mock_structure}
        
        # Act
        result = af2_predict(sequence)
        
        # Assert
        assert isinstance(result, dict)
        assert "seq_0" in result
        mock_predict.assert_called_once_with([sequence], None)

    @patch('enzy_htp._interface.alphafold_interface.interface.alphafold.predict')
    def test_af2_predict_multiple_sequences(self, mock_predict):
        """Test af2_predict with multiple sequences."""
        # Arrange
        sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY"]
        mock_structures = {
            "seq_0": Mock(spec=Structure),
            "seq_1": Mock(spec=Structure)
        }
        mock_predict.return_value = mock_structures
        
        # Act
        result = af2_predict(sequences)
        
        # Assert
        assert isinstance(result, dict)
        assert len(result) == 2
        mock_predict.assert_called_once_with(sequences, None)

    @patch('enzy_htp._interface.alphafold_interface.interface.alphafold.predict')
    def test_af2_predict_with_cluster_config(self, mock_predict):
        """Test af2_predict with cluster job configuration."""
        # Arrange
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        cluster_config = ClusterJobConfig(res_keywords={"partition": "gpu"})
        mock_structure = Mock(spec=Structure)
        mock_predict.return_value = {"seq_0": mock_structure}
        
        # Act
        result = af2_predict(sequence, cluster_job_config=cluster_config)
        
        # Assert
        assert isinstance(result, dict)
        mock_predict.assert_called_once_with([sequence], cluster_config)
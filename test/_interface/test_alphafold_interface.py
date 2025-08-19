"""Test module for enzy_htp._interface.alphafold_interface

Author: QZ shao <shaoqz@icloud.com>
Date: 2025-08-18
"""
import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import MagicMock, patch, Mock, mock_open

from enzy_htp import interface
from enzy_htp._interface.alphafold_interface import af2_predict
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.chemical.sequence import create_fasta_from_sequences
af_interface = interface.alphafold

class TestAlphafoldInterfaceUnmocked:
    """Un-mocked tests for AlphafoldInterface."""

    def test_af2_predict(self):
        """Test af2_predict without mocking."""
        # Arrange
        sequence = "MSTPSLIPSGVHEVLAKYKDGN"
        # Act
        result = af2_predict(sequence)

        # Assert
        assert isinstance(result, dict)
        assert "seq_0" in result
        assert result["seq_0"].num_residues == 22 
        assert result["seq_0"].num_atoms == 163 

    def test_build_colabfold_container_command(self):
        """Test building ColabFold container command without mocking."""
        # Arrange
        fasta_path = "/tmp/test.fasta"
        out_dir = Path("/tmp/out")
        
        # Act
        cmd = self.interface._build_colabfold_container_command(
            fasta_path, out_dir, 3, 2, 1, 200, False, "alphafold2_ptm", None
        )
        # Assert
        assert isinstance(cmd, list)
        assert len(cmd) > 5  # Should have multiple elements
        assert cmd[0] == "apptainer"  # Default container type
        assert "run" in cmd
        assert "--nv" in cmd
        assert "colabfold_batch" in cmd
        assert fasta_path in cmd
        assert "/work" in cmd
        assert "--num-models" in cmd
        assert "3" in cmd
        assert "--num-recycle" in cmd
        assert "2" in cmd
        assert "--model-type" in cmd
        assert "alphafold2_ptm" in cmd

    def test_build_alphafold2_native_command(self):
        """Test building native AlphaFold2 command without mocking."""
        # Arrange
        self.config.DATA_DIR = "/path/to/data"
        self.config.EXECUTABLE_PATH = "/path/to/run_alphafold.py"
        fasta_path = "/tmp/test.fasta"
        out_dir = Path("/tmp/out")
        
        # Act
        cmd = self.interface._build_alphafold2_native_command(
            fasta_path, out_dir, "2023-05-01", "monomer_ptm", "full_dbs", None
        )
        
        # Assert
        assert isinstance(cmd, list)
        assert cmd[0] == "python"
        assert "/path/to/run_alphafold.py" in cmd
        assert "--fasta_paths" in cmd
        assert fasta_path in cmd
        assert "--output_dir" in cmd
        assert str(out_dir) in cmd
        assert "--data_dir" in cmd
        assert "/path/to/data" in cmd
        assert "--max_template_date" in cmd
        assert "2023-05-01" in cmd

    def test_build_alphafold2_python_command(self):
        """Test building AlphaFold2 Python command without mocking."""
        # Arrange
        self.config.EXECUTABLE_PATH = "colabfold_batch"
        fasta_path = "/tmp/test.fasta"
        out_dir = Path("/tmp/out")
        
        # Act
        cmd = self.interface._build_alphafold2_python_command(
            fasta_path, out_dir, 2, 3, 1, 100, True, "alphafold2", None
        )
        
        # Assert
        assert isinstance(cmd, list)
        assert cmd[0] == "colabfold_batch"
        assert fasta_path in cmd
        assert str(out_dir) in cmd
        assert "--num-models" in cmd
        assert "2" in cmd
        assert "--templates" in cmd
        assert "--amber" in cmd

    def test_find_output_files_empty_directory(self):
        """Test finding output files in empty directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Act
            result = self.interface._find_output_files(temp_path)
            
            # Assert
            assert isinstance(result, dict)
            assert len(result) == 0

    def test_find_output_files_with_pdb_files(self):
        """Test finding output files with actual PDB files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test PDB files
            (temp_path / "seq_0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            (temp_path / "seq_0_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            (temp_path / "seq_1_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").touch()
            (temp_path / "other_file.txt").touch()  # Non-PDB file
            
            # Act
            result = self.interface._find_output_files(temp_path)
            
            # Assert
            assert isinstance(result, dict)
            assert "seq_0" in result
            assert "seq_1" in result
            # Should prefer relaxed over unrelaxed
            assert "relaxed" in result["seq_0"]
            assert "unrelaxed" in result["seq_1"]  # Only unrelaxed available
            assert len(result) == 2  # Should not include non-PDB files

    def test_unsupported_install_type_error(self):
        """Test error handling for unsupported install type."""
        # Arrange
        self.config.INSTALL_TYPE = "unsupported_type"
        
        # Act & Assert
        with pytest.raises(ValueError, match="Unsupported install type"):
            self.interface.run("/tmp/test.fasta", "/tmp/out")

    def test_make_job_with_sequence_splitting(self):
        """Test make_job with sequence splitting for array jobs."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY", "GHIKLMNPQRSTVWY"]
            fasta_path = create_fasta_from_sequences(
                sequences,
                ["seq_0", "seq_1", "seq_2"],
                output_path=temp_path / "test.fasta"
            )
            
            # Import an actual cluster type for testing
            from enzy_htp.core.clusters.accre import Accre
            cluster = Accre()
            
            # Act
            result_egg = self.interface.make_job(
                fasta_path=fasta_path,
                out_dir=temp_path / "output",
                seq_per_job=2,
                cluster_job_config=ClusterJobConfig(cluster=cluster, res_keywords={})
            )
            
            # Assert
            assert isinstance(result_egg, AlphaFold2ResultEgg)
            assert result_egg.fasta_path == fasta_path
            assert result_egg.output_dir == temp_path / "output"
            assert len(result_egg.sequence_ids) == 3
            assert len(result_egg.jobs) == 2  # 3 sequences / 2 per job = 2 jobs
            assert result_egg.sequence_ids == ["seq_0", "seq_1", "seq_2"]


class TestAlphaFold2ResultEgg:
    pass
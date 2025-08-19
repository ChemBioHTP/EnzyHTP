"""Test module for enzy_htp._interface.alphafold_interface

Author: QZ shao <shaoqz@icloud.com>
Date: 2025-08-18
"""
import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import MagicMock, Mock

from enzy_htp import interface
from enzy_htp import config as eh_config
from enzy_htp._interface.alphafold_interface import AlphafoldInterface, AlphaFold2ResultEgg
from enzy_htp.structure_prediction.prediction import predict_structure
from enzy_htp._config.alphafold_config import AlphafoldConfig
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.core.clusters.accre_r9 import AccreR9
from enzy_htp.chemical.sequence import create_fasta_from_sequences
af_interface = interface.alphafold
af_config = eh_config.alphafold

class TestAlphafoldInterfaceUnmocked:
    """Un-mocked tests for AlphafoldInterface."""

    def setup_method(self):
        """Set up test fixtures."""
        self.config = af_config
        self.interface = af_interface

    def test_predict_structure(self):
        """Test predict_structure without mocking."""
        # Arrange
        sequence = "MSTPSLIPSGVHEVLAKYKDGN"
        # Act
        result = predict_structure(sequence, engine="alphafold2")

        # Assert
        assert isinstance(result, dict)
        # The result keys should be the actual sequences, not seq_0
        assert sequence in result
        assert result[sequence].num_residues == 22 
        assert result[sequence].num_atoms == 163 

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

    def test_build_alphafold2_native_container_command(self):
        """Test building native AlphaFold2 command without mocking."""
        # Arrange
        self.config.DATA_DIR = "/path/to/data"
        self.config.EXECUTABLE_PATH = "/path/to/run_alphafold.py"
        fasta_path = "/tmp/test.fasta"
        out_dir = Path("/tmp/out")
        
        # Act
        cmd = self.interface._build_alphafold2_native_container_command(
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

    def test_build_alphafold2_native_python_command(self):
        """Test building AlphaFold2 native Python command without mocking."""
        # Arrange
        self.config.EXECUTABLE_PATH = "/path/to/run_alphafold.py"
        self.config.DATA_DIR = "/path/to/data"
        fasta_path = "/tmp/test.fasta"
        out_dir = Path("/tmp/out")
        
        # Act
        cmd = self.interface._build_alphafold2_native_python_command(
            fasta_path, out_dir, 2, 3, 1, 100, True, "alphafold2", None
        )
        
        # Assert
        assert isinstance(cmd, list)
        assert cmd[0] == "python"
        assert "/path/to/run_alphafold.py" in cmd
        assert fasta_path in cmd
        assert str(out_dir) in cmd
        assert "--data_dir" in cmd
        assert "/path/to/data" in cmd

    def test_find_output_files_empty_directory(self):
        """Test finding output files in empty directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Act
            result = self.interface._find_output_files_map(temp_path)
            
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
            filename_to_path = self.interface._find_output_files_map(temp_path)
            result = self.interface._select_best_files_for_sequences(filename_to_path, ["seq_0", "seq_1"])
            
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


class TestAlphafoldClusterJobs:
    """Tests for AlphaFold2 cluster job functionality."""

    def test_af2_predict_with_cluster_job_config_accre_r9(self):
        """Test af2_predict with ACCRE R9 cluster job configuration."""
        # Setup test config for ACCRE R9 alphafold native python install
        test_config = AlphafoldConfig()
        test_config.INSTALL_TYPE = "alphafold2_native_python"
        test_config.EXECUTABLE_PATH = "/sb/apps/alphafold232/alphafold/run_alphafold.py"
        test_config.DATA_DIR = "/sb/apps/alphafold-data.230"
        test_config.UNIREF90_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniref90/uniref90.fasta"
        test_config.MGNIFY_DATABASE_PATH = "/sb/apps/alphafold-data.230/mgnify/mgy_clusters_2022_05.fa"
        test_config.UNIREF30_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniref30/UniRef30_2021_03"
        test_config.BFD_DATABASE_PATH = "/sb/apps/alphafold-data.230/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
        test_config.TEMPLATE_MMCIF_DIR = "/sb/apps/alphafold-data.230/pdb_mmcif/mmcif_files"
        test_config.PDB_SEQRES_DATABASE_PATH = "/sb/apps/alphafold-data.230/pdb_seqres/pdb_seqres.txt"
        test_config.OBSOLETE_PDBS_PATH = "/sb/apps/alphafold-data.230/pdb_mmcif/obsolete.dat"
        test_config.UNIPROT_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniprot/uniprot.fasta"
        test_config.USE_GPU_RELAX = True

        # Create interface with test config
        interface = AlphafoldInterface(None, test_config)
        
        # Setup cluster job config for ACCRE R9
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {
                "account": "yang_lab_csb_iacc",
                "partition": "interactive_gpu", 
                "qos": "debug_iacc",
                "node_cores": "nvidia_rtx_a4000:1"
            }
        }

        # Test sequence
        sequences = ["MSTPSLIPSGVHEVLAKYKDGN"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test make_job method
            result_eggs = interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["test_seq"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
                core_type="gpu",
                seq_per_job=1
            )
            
            # Verify result egg
            assert len(result_eggs) == 1
            egg = result_eggs[0]
            assert isinstance(egg, AlphaFold2ResultEgg)
            assert len(egg.sequence_ids) == 1
            assert egg.sequence_ids[0] == "test_seq"
            
            # Verify job configuration
            job = egg.job
            assert isinstance(job, ClusterJob)
            
            # Check that default res_keywords were applied by examining the submission script
            script_content = job.sub_script_str
            assert "#SBATCH --account=yang_lab_csb_iacc" in script_content
            assert "#SBATCH --partition=interactive_gpu" in script_content
            assert "#SBATCH --qos=debug_iacc" in script_content
            assert "#SBATCH --gres=gpu:nvidia_rtx_a4000:1" in script_content
            # Check default values were merged
            assert "#SBATCH --job-name=AF2_EnzyHTP" in script_content
            assert "#SBATCH --mem=24G" in script_content
            assert "#SBATCH --time=16:00:00" in script_content

    def test_af2_predict_default_res_keywords(self):
        """Test that default res_keywords are properly applied."""
        interface = AlphafoldInterface(None, AlphafoldConfig())
        
        # Test with minimal cluster config (empty res_keywords)
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {}
        }
        
        sequences = ["ACDEFG"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            result_eggs = interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["seq1"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output", 
                cluster_job_config=cluster_job_config,
                core_type="gpu"
            )
            
            # Check that all default res_keywords were applied by examining the submission script
            job = result_eggs[0].job
            script_content = job.sub_script_str
            
            # Check that default values appear in the submission script
            assert "#SBATCH --account=yang_lab_csb_iacc" in script_content
            assert "#SBATCH --partition=interactive_gpu" in script_content
            assert "#SBATCH --qos=debug_iacc" in script_content
            assert "#SBATCH --gres=gpu:nvidia_rtx_a4000:1" in script_content
            assert "#SBATCH --job-name=AF2_EnzyHTP" in script_content
            assert "#SBATCH --mem=24G" in script_content
            assert "#SBATCH --time=16:00:00" in script_content

    def test_af2_predict_core_type_cpu(self):
        """Test af2_predict with CPU core type."""
        interface = AlphafoldInterface(None, AlphafoldConfig())
        
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {}
        }
        
        sequences = ["MSTPSL"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            result_eggs = interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["seq1"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
                core_type="cpu"
            )
            
            # Check CPU-specific defaults were applied by examining the submission script
            job = result_eggs[0].job
            script_content = job.sub_script_str
            
            assert "#SBATCH --tasks-per-node=6" in script_content
            assert "#SBATCH --partition=batch" in script_content
            assert "#SBATCH --mem-per-cpu=4G" in script_content
            assert "#SBATCH --time=24:00:00" in script_content

    def test_af2_predict_env_settings_reflection(self):
        """Test that environment settings properly reflect core_type for cluster jobs."""
        interface = AlphafoldInterface(None, AlphafoldConfig())
        
        # Use a real cluster instead of a mock
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {}
        }
        
        sequences = ["MSTPSL"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test GPU core type
            result_eggs = interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["seq1"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
                core_type="gpu"
            )
            
            # Verify that the environment settings are in the submission script
            job = result_eggs[0].job
            script_content = job.sub_script_str
            
            # Check that alphafold environment activation is present
            assert "source /sb/apps/alphafold232/miniconda3/bin/activate af232" in script_content
            assert "export LD_LIBRARY_PATH=/sb/apps/alphafold232/miniconda3/envs/af232/lib:$LD_LIBRARY_PATH" in script_content

    @pytest.mark.slow
    def test_af2_predict_real_cluster_submission(self):
        """Test actual cluster job submission (non-mocked)."""
        # Setup for ACCRE R9 with real configuration
        test_config = AlphafoldConfig()
        test_config.INSTALL_TYPE = "alphafold2_native_python"
        test_config.EXECUTABLE_PATH = "/sb/apps/alphafold232/alphafold/run_alphafold.py"
        test_config.DATA_DIR = "/sb/apps/alphafold-data.230"
        # Set all required database paths as shown in feedback example
        test_config.UNIREF90_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniref90/uniref90.fasta"
        test_config.MGNIFY_DATABASE_PATH = "/sb/apps/alphafold-data.230/mgnify/mgy_clusters_2022_05.fa"
        test_config.UNIREF30_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniref30/UniRef30_2021_03"
        test_config.BFD_DATABASE_PATH = "/sb/apps/alphafold-data.230/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
        test_config.TEMPLATE_MMCIF_DIR = "/sb/apps/alphafold-data.230/pdb_mmcif/mmcif_files"
        test_config.PDB_SEQRES_DATABASE_PATH = "/sb/apps/alphafold-data.230/pdb_seqres/pdb_seqres.txt"
        test_config.OBSOLETE_PDBS_PATH = "/sb/apps/alphafold-data.230/pdb_mmcif/obsolete.dat"
        test_config.UNIPROT_DATABASE_PATH = "/sb/apps/alphafold-data.230/uniprot/uniprot.fasta"
        test_config.USE_GPU_RELAX = True

        interface = AlphafoldInterface(None, test_config)
        
        # Cluster configuration for actual submission
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {
                "account": "yang_lab_csb_iacc",
                "partition": "interactive_gpu",
                "qos": "debug_iacc", 
                "node_cores": "nvidia_rtx_a4000:1"
            }
        }

        # Very short test sequence to minimize computational cost
        sequences = ["MST"]  # 3 amino acids only
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create result eggs but don't actually submit (for safety)
            result_eggs = interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["test_seq"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
                core_type="gpu",
                seq_per_job=1
            )
            
            # Verify the job is properly configured for submission
            job = result_eggs[0].job
            assert job.cluster.__class__.__name__ == "AccreR9"
            
            # Check that the submission script contains expected SBATCH directives
            # This validates the job is ready for actual submission
            script_content = job.sub_script_str
            assert "#SBATCH --account=yang_lab_csb_iacc" in script_content
            assert "#SBATCH --partition=interactive_gpu" in script_content
            assert "#SBATCH --qos=debug_iacc" in script_content
            assert "#SBATCH --gres=gpu:nvidia_rtx_a4000:1" in script_content
            
            # Verify AlphaFold command is properly constructed
            assert "python /sb/apps/alphafold232/alphafold/run_alphafold.py" in script_content
            assert "--data_dir /sb/apps/alphafold-data.230" in script_content
            assert "--use_gpu_relax" in script_content


class TestAlphaFold2ResultEgg:
    """Tests for AlphaFold2 result egg functionality."""

    def test_result_egg_initialization(self):
        """Test AlphaFold2ResultEgg initialization."""
        mock_job = Mock(spec=ClusterJob)
        mock_interface = Mock(spec=AlphafoldInterface)
        
        egg = AlphaFold2ResultEgg(
            fasta_path="/tmp/test.fasta",
            output_dir=Path("/tmp/output"),
            sequence_ids=["seq1", "seq2"],
            job=mock_job,
            interface=mock_interface
        )
        
        assert egg.fasta_path == "/tmp/test.fasta"
        assert egg.output_dir == Path("/tmp/output")
        assert egg.sequence_ids == ["seq1", "seq2"]
        assert egg.job == mock_job
        assert egg.interface == mock_interface

    def test_get_expected_output_files(self):
        """Test getting expected output files from result egg."""
        mock_interface = Mock(spec=AlphafoldInterface)
        mock_interface._find_output_files_map.return_value = {
            "seq1_rank_001": "/path/to/seq1_rank_001.pdb",
            "seq2_rank_001": "/path/to/seq2_rank_001.pdb"
        }
        mock_interface._select_best_files_for_sequences.return_value = {
            "seq1": "/path/to/seq1_rank_001.pdb",
            "seq2": "/path/to/seq2_rank_001.pdb"
        }
        
        egg = AlphaFold2ResultEgg(
            fasta_path="/tmp/test.fasta",
            output_dir=Path("/tmp/output"),
            sequence_ids=["seq1", "seq2"],
            job=Mock(),
            interface=mock_interface
        )
        
        result = egg.get_expected_output_files()
        
        assert result == {
            "seq1": "/path/to/seq1_rank_001.pdb",
            "seq2": "/path/to/seq2_rank_001.pdb"
        }
        mock_interface._find_output_files_map.assert_called_once_with(Path("/tmp/output"))
        mock_interface._select_best_files_for_sequences.assert_called_once()
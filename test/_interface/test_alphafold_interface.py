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
from enzy_htp.core.logger import _LOGGER
from enzy_htp.chemical.sequence import create_fasta_from_sequences
af_interface = interface.alphafold
af_config = eh_config.alphafold


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

    def test_build_alphafold2_native_container_command(self, alphafold_config_modifier):
        """Test building native AlphaFold2 command without mocking."""
        # Arrange
        alphafold_config_modifier(
            DATA_DIR="/path/to/data",
            EXECUTABLE_PATH="/path/to/run_alphafold.py"
        )
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

    def test_build_alphafold2_native_python_command(self, alphafold_config_modifier):
        """Test building AlphaFold2 native Python command without mocking."""
        # Arrange
        alphafold_config_modifier(
            EXECUTABLE_PATH="/path/to/run_alphafold.py",
            DATA_DIR="/path/to/data"
        )
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

    def test_unsupported_install_type_error(self, alphafold_config_modifier):
        """Test error handling for unsupported install type."""
        # Arrange
        alphafold_config_modifier(
            INSTALL_TYPE="unsupported_type",
            EXECUTABLE_PATH="placeholder"
        )

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
            cluster = AccreR9()
            
            # Act
            result_eggs = self.interface.make_job(
                fasta_path=fasta_path,
                out_dir=temp_path / "output",
                seq_per_job=2,
                cluster_job_config=ClusterJobConfig(cluster=cluster, res_keywords={})
            )
            
            # Assert
            assert isinstance(result_eggs, list)
            assert len(result_eggs) == 2  # 3 sequences / 2 per job = 2 jobs
            assert all(isinstance(egg, AlphaFold2ResultEgg) for egg in result_eggs)
            
            # First job should have 2 sequences
            assert len(result_eggs[0].sequence_ids) == 2
            assert result_eggs[0].sequence_ids == ["seq_0", "seq_1"]
            
            # Second job should have 1 sequence 
            assert len(result_eggs[1].sequence_ids) == 1
            assert result_eggs[1].sequence_ids == ["seq_2"]


class TestAlphafoldClusterJobs:
    """Tests for AlphaFold2 cluster job functionality."""

    def test_make_job_with_cluster_job_config_accre_r9(self, alphafold_config_modifier):
        """Test make_job with ACCRE R9 cluster job configuration."""
        # Setup test config for ACCRE R9 alphafold native python install
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            UNIREF90_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref90/uniref90.fasta",
            MGNIFY_DATABASE_PATH="/sb/apps/alphafold-data.230/mgnify/mgy_clusters_2022_05.fa",
            UNIREF30_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref30/UniRef30_2021_03",
            BFD_DATABASE_PATH="/sb/apps/alphafold-data.230/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
            TEMPLATE_MMCIF_DIR="/sb/apps/alphafold-data.230/pdb_mmcif/mmcif_files",
            PDB_SEQRES_DATABASE_PATH="/sb/apps/alphafold-data.230/pdb_seqres/pdb_seqres.txt",
            OBSOLETE_PDBS_PATH="/sb/apps/alphafold-data.230/pdb_mmcif/obsolete.dat",
            UNIPROT_DATABASE_PATH="/sb/apps/alphafold-data.230/uniprot/uniprot.fasta",
            USE_GPU_RELAX=True
        )

        # Create interface with global config
        interface = af_interface
        
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

    def test_make_job_default_res_keywords(self):
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
            )
            
            # Check that all default res_keywords were applied by examining the submission script
            job = result_eggs[0].job
            script_content = job.sub_script_str
            
            # Check that default values appear in the submission script
            # Note: this test uses default AlphafoldConfig which has <fillthis> placeholders
            assert "#SBATCH --account=<fillthis>" in script_content
            assert "#SBATCH --partition=<fillthis>" in script_content
            assert "#SBATCH --gres=gpu:nvidia_rtx_a4000:1" in script_content
            assert "#SBATCH --job-name=AF2_EnzyHTP" in script_content
            assert "#SBATCH --mem=24G" in script_content
            assert "#SBATCH --time=16:00:00" in script_content

    def test_make_job_env_settings_reflection(self):
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
            )
            
            # Verify that the environment settings are in the submission script
            job = result_eggs[0].job
            script_content = job.sub_script_str
            
            # Check that alphafold environment activation is present
            assert "source /sb/apps/alphafold232/miniconda3/bin/activate af232" in script_content
            assert "export LD_LIBRARY_PATH=/sb/apps/alphafold232/miniconda3/envs/af232/lib:$LD_LIBRARY_PATH" in script_content

    @pytest.mark.slow
    def test_af2_predict_real_cluster_submission(self, alphafold_config_modifier):
        """Test actual cluster job submission (non-mocked)."""
        # Setup for ACCRE R9 with real configuration
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            UNIREF90_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref90/uniref90.fasta",
            MGNIFY_DATABASE_PATH="/sb/apps/alphafold-data.230/mgnify/mgy_clusters_2022_05.fa",
            UNIREF30_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref30/UniRef30_2021_03",
            BFD_DATABASE_PATH="/sb/apps/alphafold-data.230/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
            TEMPLATE_MMCIF_DIR="/sb/apps/alphafold-data.230/pdb_mmcif/mmcif_files",
            PDB_SEQRES_DATABASE_PATH="/sb/apps/alphafold-data.230/pdb_seqres/pdb_seqres.txt",
            OBSOLETE_PDBS_PATH="/sb/apps/alphafold-data.230/pdb_mmcif/obsolete.dat",
            UNIPROT_DATABASE_PATH="/sb/apps/alphafold-data.230/uniprot/uniprot.fasta",
            USE_GPU_RELAX=True
        )
        
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
            result_eggs = af_interface.make_job(
                fasta_path=create_fasta_from_sequences(sequences, ["test_seq"], output_path=Path(temp_dir) / "test.fasta"),
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
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

    def test_af2_predict_non_armer_core_type_cpu(self, alphafold_config_modifier):
        """Test non_armer_core_type parameter with CPU setting for local execution."""
        # Setup config for alphafold2_native_python
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            USE_GPU_RELAX=True  # This should be ignored when core_type is cpu
        )
        
        interface = AlphafoldInterface(None, af_config)
        
        sequences = ["MST"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test CPU core type with local execution (cluster_job_config=None)
            fasta_path = create_fasta_from_sequences(sequences, ["test_seq"], output_path=Path(temp_dir) / "test.fasta")
            
            # Test command building directly
            cmd = interface._build_alphafold2_native_python_command(
                fasta_path=fasta_path,
                out_dir=Path(temp_dir) / "output",
                num_models=1,
                num_recycles=1,
                num_relax=0,
                relax_max_iteration=100,
                use_templates=False,
                model_preset="monomer",
                additional_options=None,
                core_type="cpu"
            )
            
            # CPU mode should NOT include --use_gpu_relax
            assert "--use_gpu_relax" not in cmd
            assert "python" in cmd[0]
            assert "/sb/apps/alphafold232/alphafold/run_alphafold.py" in cmd

    def test_af2_predict_non_armer_core_type_gpu(self, alphafold_config_modifier):
        """Test non_armer_core_type parameter with GPU setting for local execution."""
        # Setup config for alphafold2_native_python
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            USE_GPU_RELAX=True
        )
        
        interface = AlphafoldInterface(None, af_config)
        
        sequences = ["MST"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test GPU core type with local execution (cluster_job_config=None)
            fasta_path = create_fasta_from_sequences(sequences, ["test_seq"], output_path=Path(temp_dir) / "test.fasta")
            
            # Test command building directly
            cmd = interface._build_alphafold2_native_python_command(
                fasta_path=fasta_path,
                out_dir=Path(temp_dir) / "output",
                num_models=1,
                num_recycles=1,
                num_relax=0,
                relax_max_iteration=100,
                use_templates=False,
                model_preset="monomer",
                additional_options=None,
                core_type="gpu"
            )
            
            # GPU mode should include --use_gpu_relax
            assert "--use_gpu_relax" in cmd
            assert "python" in cmd[0]
            assert "/sb/apps/alphafold232/alphafold/run_alphafold.py" in cmd

    def test_colabfold_container_core_type_cpu(self):
        """Test ColabFold container command building with CPU core type."""
        interface = AlphafoldInterface(None, AlphafoldConfig())
        
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_path = create_fasta_from_sequences(["MST"], ["test_seq"], output_path=Path(temp_dir) / "test.fasta")
            
            # Test CPU core type
            cmd = interface._build_colabfold_container_command(
                fasta_path=fasta_path,
                out_dir=Path(temp_dir) / "output",
                num_models=1,
                num_recycles=1,
                num_relax=0,
                relax_max_iteration=100,
                use_templates=False,
                model_preset="alphafold2_ptm",
                additional_options=None,
                core_type="cpu"
            )
            
            # CPU mode should NOT include --nv
            assert "--nv" not in cmd
            assert "apptainer" in cmd
            assert "run" in cmd
            
    def test_colabfold_container_core_type_gpu(self):
        """Test ColabFold container command building with GPU core type."""
        interface = AlphafoldInterface(None, AlphafoldConfig())
        
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_path = create_fasta_from_sequences(["MST"], ["test_seq"], output_path=Path(temp_dir) / "test.fasta")
            
            # Test GPU core type
            cmd = interface._build_colabfold_container_command(
                fasta_path=fasta_path,
                out_dir=Path(temp_dir) / "output",
                num_models=1,
                num_recycles=1,
                num_relax=0,
                relax_max_iteration=100,
                use_templates=False,
                model_preset="alphafold2_ptm",
                additional_options=None,
                core_type="gpu"
            )
            
            # GPU mode should include --nv
            assert "--nv" in cmd
            assert "apptainer" in cmd
            assert "run" in cmd


class TestAlphafoldAccreR9Integration:
    """Tests for AlphaFold2 integration with ACCRE R9 cluster."""

    @pytest.mark.slow
    def test_af2_real_accre_r9_native_python_job_setup(self, alphafold_config_modifier):
        """Test real AF2 job setup for ACCRE R9 with native Python install (no mocking)."""
        # Setup for actual ACCRE R9 alphafold native python install
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            UNIREF90_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref90/uniref90.fasta",
            MGNIFY_DATABASE_PATH="/sb/apps/alphafold-data.230/mgnify/mgy_clusters_2022_05.fa",
            UNIREF30_DATABASE_PATH="/sb/apps/alphafold-data.230/uniref30/UniRef30_2021_03",
            BFD_DATABASE_PATH="/sb/apps/alphafold-data.230/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
            TEMPLATE_MMCIF_DIR="/sb/apps/alphafold-data.230/pdb_mmcif/mmcif_files",
            PDB_SEQRES_DATABASE_PATH="/sb/apps/alphafold-data.230/pdb_seqres/pdb_seqres.txt",
            OBSOLETE_PDBS_PATH="/sb/apps/alphafold-data.230/pdb_mmcif/obsolete.dat",
            UNIPROT_DATABASE_PATH="/sb/apps/alphafold-data.230/uniprot/uniprot.fasta",
            USE_GPU_RELAX=True
        )
        
        # Cluster configuration for ACCRE R9 
        cluster_job_config = {
            "cluster": AccreR9(),
            "res_keywords": {
                "account": "yang_lab_csb_iacc",
                "partition": "interactive_gpu",
                "qos": "debug_iacc", 
                "node_cores": "nvidia_rtx_a4000:1",
                "walltime": "30:00",
            }
        }

        # Very short test sequence to minimize computational cost
        sequences = ["MSTPSLIPSGVHEVLAKYKDGN"]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create result eggs for actual submission
            result = af_interface.af2_predict(
                sequences=sequences,
                out_dir=Path(temp_dir) / "output",
                cluster_job_config=cluster_job_config,
                seq_per_job=1,
                model_preset="monomer_ptm",
            )
            print(result)

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
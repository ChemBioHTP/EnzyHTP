"""Test module for enzy_htp._interface.alphafold_interface

Author: QZ shao <shaoqz@icloud.com>
Date: 2025-08-18
"""
import pytest
import tempfile
import os
from pathlib import Path

from enzy_htp import interface
from enzy_htp import config as eh_config
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.core.clusters.accre_r9 import AccreR9
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


class TestAlphafoldInterface:
    """Un-mocked tests for AlphafoldInterface."""

    def setup_method(self):
        """Set up test fixtures."""
        self.config = af_config
        self.interface = af_interface

    def test_af2_predict_colabfold(self):
        """Test predict_structure without mocking use colabfold"""
        test_output_dir = Path(f"{WORK_DIR}/test_af2_output")
        # Arrange
        sequence = "MSTPSLIPSGVHEVLAKYKDGN"
        # Act
        result = self.interface.af2_predict(sequence, work_dir=test_output_dir)

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

    def test_af2_predict_colabfold_multimer(self):
        """Test predict_structure without mocking use colabfold"""
        test_output_dir = Path(f"{WORK_DIR}/test_af2_output")
        # Arrange
        sequence = [["MSTPSLIPSGVHEVLAKYKDGN", "MSTPSLIPSGVHEVLAKYKDGN"],]
        # Act
        result = self.interface.af2_predict(sequence, work_dir=test_output_dir)

        # Assert
        assert isinstance(result, dict)
        # The result keys should be the actual sequences, not seq_0
        assert tuple(sequence[0]) in result
        
        # Check comprehensive output structure
        seq_result = result[tuple(sequence[0])]
        assert isinstance(seq_result, dict)
        
        # Check best model
        assert "best_model" in seq_result
        assert "best_model_plddt" in seq_result
        assert "best_model_index" in seq_result
        assert seq_result["best_model"].num_residues == 44
        assert seq_result["best_model"].num_atoms == 326
        
        # Check individual models exist
        assert "model_1" in seq_result or "model_2" in seq_result or "model_3" in seq_result or "model_4" in seq_result or "model_5" in seq_result
        
        # Check pLDDT scores are valid
        assert isinstance(seq_result["best_model_plddt"], list)
        assert len(seq_result["best_model_plddt"]) == 44  # Same as number of residues
        fs.safe_rmdir(test_output_dir)

    def test_parse_native_af2_results(self):
        """Test parsing native AlphaFold2 results format using reference data."""
        # Arrange
        ref_output_dir = Path(f"{DATA_DIR}/test_af2_output_ref")
        seq_id = "seq_0"
        parser = PDBParser()
        
        # Act
        result = self.interface._parse_comprehensive_results(ref_output_dir, seq_id, parser)
        
        # Assert
        assert isinstance(result, dict)
        assert len(result) > 0, "Should have parsed some results"
        
        # Check comprehensive output structure
        assert "best_model" in result
        assert "best_model_plddt" in result
        assert "best_model_index" in result
        
        # Verify best model structure
        assert result["best_model"].num_residues == 22
        assert result["best_model"].num_atoms == 163
        
        # Check pLDDT scores are valid
        assert isinstance(result["best_model_plddt"], list)
        assert len(result["best_model_plddt"]) == 22  # Same as number of residues
        
        # Based on ranking_debug.json, model_5 should be the best
        assert result["best_model_index"] == 5
        
        # Check that all 5 models are present
        for i in range(1, 6):
            assert f"model_{i}" in result
            assert f"model_{i}_plddt" in result
            assert isinstance(result[f"model_{i}_plddt"], list)
            assert len(result[f"model_{i}_plddt"]) == 22
        
        # Verify no PAE in output (removed as per feedback)
        assert "best_model_pae" not in result
        for i in range(1, 6):
            assert f"model_{i}_pae" not in result

    def test_parse_native_af2_multimer_results(self):
        """Test parsing native AlphaFold2 results format using reference data."""
        # Arrange
        ref_output_dir = Path(f"{DATA_DIR}/test_af2_output_multimer_ref")
        seq_id = "seq_0"
        parser = PDBParser()
        
        # Act
        result = self.interface._parse_comprehensive_results(ref_output_dir, seq_id, parser)
        
        # Assert
        assert isinstance(result, dict)
        assert len(result) > 0, "Should have parsed some results"
        
        # Check comprehensive output structure
        assert "best_model" in result
        assert "best_model_plddt" in result
        assert "best_model_index" in result
        
        # Verify best model structure
        assert result["best_model"].num_residues == 44
        assert result["best_model"].num_atoms == 326
        
        # Check pLDDT scores are valid
        assert isinstance(result["best_model_plddt"], list)
        assert len(result["best_model_plddt"]) == 44  # Same as number of residues
        
        # Based on ranking_debug.json, model_5 should be the best
        assert result["best_model_index"] == 2
        
        # Check that all 5 models are present
        for i in range(1, 6):
            assert f"model_{i}" in result
            assert f"model_{i}_plddt" in result
            assert isinstance(result[f"model_{i}_plddt"], list)
            assert len(result[f"model_{i}_plddt"]) == 44
        
        # Verify no PAE in output (removed as per feedback)
        assert "best_model_pae" not in result
        for i in range(1, 6):
            assert f"model_{i}_pae" not in result

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

class TestAlphafoldInterfaceOnAccreR9:
    """Tests for AlphaFold2 integration with ACCRE R9 cluster."""

    @pytest.mark.slow
    def test_af2_real_accre_r9_native_python_job_setup(self, alphafold_config_modifier):
        """Test real AF2 job setup for ACCRE R9 with native Python install (no mocking).
        The test take 13~27 mins even with precomputed msas"""
        # Setup for actual ACCRE R9 alphafold native python install
        dummydb_path = f"{DATA_DIR}/dummy_af2_template_database"
        dummy_msa_path = f"{DATA_DIR}/dummy_msas_1"
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            PDB_SEQRES_DATABASE_PATH=f"{dummydb_path}/dummy_fas.fas",
            TEMPLATE_MMCIF_DIR=f"{dummydb_path}/",
            PDB70_DATABASE_PATH=f"{dummydb_path}/dummydb",
            OBSOLETE_PDBS_PATH=f"{dummydb_path}/dummy_obsolete.dat",
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
        test_output_dir = Path(f"{WORK_DIR}/test_af2_output")
        fs.safe_mkdir(str(test_output_dir))
        fs.safe_mkdir(str(test_output_dir / "seq_0"))
        fs.safe_cp(dummy_msa_path, f"{test_output_dir}/seq_0/msas", allow_dir=True)

        try:
            # Create result eggs for actual submission
            result = af_interface.af2_predict(
                sequences=sequences,
                work_dir=test_output_dir,
                cluster_job_config=cluster_job_config,
                seq_per_job=1,
                model_preset="monomer_ptm",
                use_precomputed_msas=True,
            )
            assert result
            assert isinstance(result, dict)
            for seq in sequences:
                assert seq in result
                seq_result = result[seq]
                assert isinstance(seq_result, dict)
                assert "best_model" in seq_result
                assert "best_model_plddt" in seq_result
                assert "best_model_index" in seq_result
                assert isinstance(seq_result["best_model_plddt"], list)
                assert len(seq_result["best_model_plddt"]) == 22
                assert isinstance(seq_result["best_model"], Structure)
        finally:
            # Cleanup test directory
            if test_output_dir.exists():
                fs.safe_rmdir(test_output_dir)

    @pytest.mark.slow
    def test_af2_real_accre_r9_native_python_multi_job(self, alphafold_config_modifier):
        """Test real AF2 job setup for ACCRE R9 with native Python install (no mocking).
        The test take 13~27 mins even with precomputed msas"""
        # Setup for actual ACCRE R9 alphafold native python install
        dummydb_path = f"{DATA_DIR}/dummy_af2_template_database"
        dummy_msa_path = f"{DATA_DIR}/dummy_msas_1"
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            PDB_SEQRES_DATABASE_PATH=f"{dummydb_path}/dummy_fas.fas",
            TEMPLATE_MMCIF_DIR=f"{dummydb_path}/",
            PDB70_DATABASE_PATH=f"{dummydb_path}/dummydb",
            OBSOLETE_PDBS_PATH=f"{dummydb_path}/dummy_obsolete.dat",
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
        sequences = ["MSTPSLIPSGVHEVLAKYKDGN", "MSTPSLAAAGVHEVLAKYKDGN"]
        test_output_dir = Path(f"{WORK_DIR}/test_af2_output")
        fs.safe_mkdir(str(test_output_dir))
        fs.safe_mkdir(str(test_output_dir / "seq_0"))
        fs.safe_mkdir(str(test_output_dir / "seq_1"))
        fs.safe_cpdir(dummy_msa_path, f"{test_output_dir}/seq_0/msas/")
        fs.safe_cpdir(dummy_msa_path, f"{test_output_dir}/seq_1/msas/")

        try:
            # Create result eggs for actual submission
            result = af_interface.af2_predict(
                sequences=sequences,
                work_dir=test_output_dir,
                cluster_job_config=cluster_job_config,
                seq_per_job=1,
                model_preset="monomer_ptm",
                use_precomputed_msas=True,
            )
            assert result
            assert isinstance(result, dict)
            for seq in sequences:
                assert seq in result
                seq_result = result[seq]
                assert isinstance(seq_result, dict)
                assert "best_model" in seq_result
                assert "best_model_plddt" in seq_result
                assert "best_model_index" in seq_result
                assert isinstance(seq_result["best_model_plddt"], list)
                assert len(seq_result["best_model_plddt"]) == 22
                assert isinstance(seq_result["best_model"], Structure)
        finally:
            # Cleanup test directory
            if test_output_dir.exists():
                fs.safe_rmdir(test_output_dir)

    @pytest.mark.slow
    def test_af2_real_accre_r9_native_python_multimer(self, alphafold_config_modifier):
        """Test real AF2 job setup for ACCRE R9 with native Python install (no mocking)."""
        # Setup for actual ACCRE R9 alphafold native python install
        dummydb_path = f"{DATA_DIR}/dummy_af2_template_database"
        dummy_msa_path = f"{DATA_DIR}/dummy_msas_1"
        alphafold_config_modifier(
            INSTALL_TYPE="alphafold2_native_python",
            EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py",
            DATA_DIR="/sb/apps/alphafold-data.230",
            PDB_SEQRES_DATABASE_PATH=f"{dummydb_path}/dummy_fas.fas",
            TEMPLATE_MMCIF_DIR=f"{dummydb_path}/",
            PDB70_DATABASE_PATH=f"{dummydb_path}/dummydb",
            OBSOLETE_PDBS_PATH=f"{dummydb_path}/dummy_obsolete.dat",
            USE_GPU_RELAX=True
        )
        
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
        sequences = [["MSTPSLIPSGVHEVLAKYKDGN", "MSTPSLIPSGVHEVLAKYKDGN"],]
        test_output_dir = Path(f"{WORK_DIR}/test_af2_output")
        fs.safe_mkdir(str(test_output_dir / "seq_0/msas"))
        fs.safe_cpdir(dummy_msa_path, f"{test_output_dir}/seq_0/msas/A/")
        fs.safe_cpdir(dummy_msa_path, f"{test_output_dir}/seq_0/msas/B/")

        try:
            # Create result eggs for actual submission
            result = af_interface.af2_predict(
                sequences=sequences,
                work_dir=test_output_dir,
                cluster_job_config=cluster_job_config,
                seq_per_job=1,
                use_precomputed_msas=True,
                num_multimer_predictions_per_model=1,
            )
            print(result)
        finally:
            # Cleanup test directory
            if test_output_dir.exists():
                fs.safe_rmdir(test_output_dir)

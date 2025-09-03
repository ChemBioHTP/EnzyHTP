"""Science API for structure prediction.

This module exposes a thin wrapper around engine-specific predictors
such as AlphaFold2. It normalizes various sequence input formats and
forwards engine-specific options.
"""
from __future__ import annotations
from typing import Any, Dict, Callable, Optional, Union, List, Tuple
from pathlib import Path

from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp.core.logger import _LOGGER
from enzy_htp.chemical.sequence import parse_fasta_file
from enzy_htp.structure import Structure
from enzy_htp import interface

def predict_structure(
        sequences: Union[str, List[str], List[List[str]], Path],
        engine: str = "alphafold2",
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
        **kwargs: Any,
    ) -> Dict[Union[str, Tuple[str, ...]], Dict[str, Any]]:
    """Predict protein structure(s) using the specified engine.

    This function standardizes sequence input and forwards options to the
    selected prediction engine (currently "alphafold2"). For AlphaFold2,
    it supports both monomers and multimers and returns comprehensive
    per-sequence results (best model, per-model structures/scores, etc.).

    Args:
        sequences: Input sequences to predict. Supported forms:
          - str: Single amino-acid sequence (monomer)
          - List[str]: Multiple monomer sequences [seq1, seq2, ...]
          - List[List[str]]: Multimers [[A_seq, B_seq], [A_seq, B_seq, C_seq], ...]
          - Path: Path to a FASTA file containing one or more sequences
        engine: Name of prediction engine. Keys of `PREDICTION_ENGINES`.
        cluster_job_config: Cluster submission configuration for engine runs.
        **kwargs: Engine-specific options forwarded as-is. For AlphaFold2
          (`interface.alphafold.af2_predict`) the commonly used options are:
          - work_dir (str | Path | None): Output/work directory
          - non_armer_core_type (str): "gpu" or "cpu" for local runs
          - cluster_job_config (ClusterJobConfig | dict | None): Cluster config
          - array_size (int): Max concurrent array jobs (cluster)
          - job_check_period (int): Seconds between job state checks
          - seq_per_job (int): Number of sequences per job (array)
          - model_preset (str | None): AF2 model preset
          - num_models (int): Number of models to generate
          - num_recycles (int): Recycling iterations
          - num_multimer_predictions_per_model (int): Multimer predictions/model
          - num_relax (int): Top-ranked structures to relax
          - relax_max_iteration (int): Max relaxation iterations
          - db_preset (str): Database preset (native AF2)
          - use_precomputed_msas (bool): Use precomputed MSAs (native AF2)
          - use_templates (bool): Whether to use templates
          - max_template_date (str | None): Max template date (native AF2)
          - additional_options (List[str] | None): Extra CLI flags
          - random_seed (int | None): Random seed for reproducibility

    Returns:
        Dict mapping each input (monomer str or multimer tuple[str, ...]) to a
        comprehensive result dict, typically including (for AlphaFold2):
          - "best_model": Structure of best-ranked model
          - "best_model_plddt": List[float] per-residue pLDDT of best model
          - "best_model_index": int (1-based index of best model)
          - "model_1", "model_1_plddt", ... for all generated models

    Raises:
        ValueError: If the specified engine is not supported or inputs invalid.
    """
    if engine not in PREDICTION_ENGINES:
        _LOGGER.error(f"Unsupported prediction engine: {engine}")
        raise ValueError(f"Unsupported engine: {engine}")

    # Parse input sequences
    parsed_sequences = _parse_sequences_input(sequences)
    
    return PREDICTION_ENGINES[engine](parsed_sequences, cluster_job_config=cluster_job_config, **kwargs)

def _parse_sequences_input(sequences: Union[str, List[str], List[List[str]], Path]) -> List[Union[str, List[str]]]:
    """Parse and normalize sequence input.

    Supports monomers and multimers. In FASTA or raw strings, a colon
    ":" within a sequence denotes multimer chains (e.g., "A:B:C"). When
    any multimer is detected, all entries are normalized to multimer form
    (List[List[str]]), with monomers wrapped as single-item lists.

    Args:
        sequences: One of
          - str: a single sequence (use ":" to separate chains for multimer)
          - List[str]: multiple monomer sequences, or strings with ":" for multimers
          - List[List[str]]: multimer definitions per target
          - Path: FASTA file path (.fasta/.fa/.fas); sequences may include ":"

    Returns:
        List of sequences in normalized form:
          - List[str] for pure monomers
          - List[List[str]] for multimers (or when any ":" is present)
    """
    if isinstance(sequences, (str, Path)):
        # Check if it's a file path
        path = Path(sequences)
        if path.exists() and path.suffix.lower() in ['.fasta', '.fa', '.fas']:
            # Parse FASTA file
            fasta_sequences = parse_fasta_file(str(path))
            # Normalize all to multimer form; wrap monomers as single-item lists
            sequences = []
            for _, seq in fasta_sequences:
                chains = seq.split(':') if ':' in seq else [seq]
                sequences.append(chains)
            return sequences
        elif isinstance(sequences, str):
            # Single sequence string
            if ':' in sequences:
                return [sequences.split(':')]
            return [sequences]
        else:
            # Path doesn't exist or isn't a FASTA file
            _LOGGER.error(f"File not found or not a FASTA file: {path}")
            raise ValueError(f"File not found or not a FASTA file: {path}")
    elif isinstance(sequences, list):
        # Already a list of sequences or list of lists (multimers)
        return sequences
    else:
        _LOGGER.error(f"Unsupported sequences input type: {type(sequences)}")
        raise ValueError(f"Unsupported sequences input type: {type(sequences)}")


PREDICTION_ENGINES: Dict[str, Callable] = {
    "alphafold2": interface.alphafold.af2_predict,
}

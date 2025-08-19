"""Science API for structure prediction."""
from __future__ import annotations
from typing import Dict, Callable, Optional, Union
from pathlib import Path

from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp.core.logger import _LOGGER
from enzy_htp.chemical.sequence import parse_fasta_file
from enzy_htp.structure import Structure
from enzy_htp import interface

def predict_structure(
        sequences: Union[str, list[str], Path], engine: str = "alphafold2", 
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None, **kwargs
    ) -> Dict[str, Structure]:
    """
    Predict protein structure(s) for the given amino acid sequence(s) using the specified engine.

    Args:
        sequences (str, list[str], or Path): A single amino acid sequence, a list of sequences, 
                  or a path to a FASTA file containing sequences to predict.
        engine (str): Name of the prediction engine to use. Supported values are keys of PREDICTION_ENGINES (default: "alphafold2").
        cluster_job_config (ClusterJobConfig or dict, optional): Configuration for submitting jobs to a cluster. Defaults to None.
        **kwargs: Engine-specific keyword arguments.

    Returns:
        dict[str, Structure]: A mapping from each input sequence (or its identifier) to the predicted Structure object.

    Raises:
        ValueError: If the specified engine is not supported.
    """
    if engine not in PREDICTION_ENGINES:
        _LOGGER.error(f"Unsupported prediction engine: {engine}")
        raise ValueError(f"Unsupported engine: {engine}")

    # Parse input sequences
    parsed_sequences = _parse_sequences_input(sequences)
    
    return PREDICTION_ENGINES[engine](parsed_sequences, cluster_job_config, **kwargs)

def _parse_sequences_input(sequences: Union[str, list[str], Path]) -> list[str]:
    """Parse sequence input from various formats.
    
    Args:
        sequences: Input sequences as string, list, or fasta file path
        
    Returns:
        List of amino acid sequences
    """
    if isinstance(sequences, (str, Path)):
        # Check if it's a file path
        path = Path(sequences)
        if path.exists() and path.suffix.lower() in ['.fasta', '.fa', '.fas']:
            # Parse FASTA file
            fasta_sequences = parse_fasta_file(str(path))
            return [seq for _, seq in fasta_sequences]
        elif isinstance(sequences, str):
            # Single sequence string
            return [sequences]
        else:
            # Path doesn't exist or isn't a FASTA file
            _LOGGER.error(f"File not found or not a FASTA file: {path}")
            raise ValueError(f"File not found or not a FASTA file: {path}")
    elif isinstance(sequences, list):
        # List of sequences
        return sequences
    else:
        _LOGGER.error(f"Unsupported sequences input type: {type(sequences)}")
        raise ValueError(f"Unsupported sequences input type: {type(sequences)}")


PREDICTION_ENGINES: Dict[str, Callable] = {
    "alphafold2": interface.alphafold.af2_predict,
}

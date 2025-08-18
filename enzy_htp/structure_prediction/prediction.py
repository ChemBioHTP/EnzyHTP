"""Science API for structure prediction."""
from __future__ import annotations
from typing import Dict, Callable, Optional, Union

from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Structure
from enzy_htp._interface.alphafold_interface import af2_predict


def predict_structure(
        sequences: Union[str, list[str]], engine: str = "alphafold2", 
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None, **kwargs
    ) -> Dict[str, Structure]:
    """
    Predict protein structure(s) for the given amino acid sequence(s) using the specified engine.

    Args:
        sequences (str or list[str]): A single amino acid sequence or a list of sequences to predict.
        engine (str): Name of the prediction engine to use. Supported values are keys of PREDICTION_ENGINES (default: "alphafold").
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

    return PREDICTION_ENGINES[engine](sequences, cluster_job_config, **kwargs)


PREDICTION_ENGINES: Dict[str, Callable] = {
    "alphafold2": af2_predict,
}

"""Science API for structure prediction."""
from __future__ import annotations

from enzy_htp.structure import Structure
from enzy_htp import interface

def predict_structure(sequence: str, engine: str = "alphafold", **kwargs) -> Structure:
    """Predict a structure from a sequence.

    Args:
        sequence: The amino acid sequence.
        engine: The prediction engine to use. Defaults to "alphafold".
        **kwargs: Additional keyword arguments for the engine.

    Returns:
        The predicted structure.
    """
    if engine == "alphafold":
        return interface.alphafold.run(sequence, **kwargs)
    else:
        raise ValueError(f"Unsupported engine: {engine}")
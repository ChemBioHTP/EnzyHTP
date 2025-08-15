"""Submodule contains code for calculating the Domain-Domain Interaction Index (DSI).

DSI is a measure of distance between two domains, defined as:
DSI = d(com1, com2) - (Rg1 + Rg2)
where d(com1, com2) is the distance between the centers of mass of the
two domains, and Rg1 and Rg2 are their respective radii of gyration.

+ dsi()
    Calculate the DSI for a trajectory from a StructureEnsemble

Author: QZ Shao <shaoqz@icloud.com>
Date: 2025-08-15
"""
from typing import List, Tuple, Callable, Dict
import numpy as np

from enzy_htp import interface
from enzy_htp.structure import StructureEnsemble
from enzy_htp.core import _LOGGER


def dsi(
    ensemble: StructureEnsemble, 
    domain1_residues: List[Tuple[str, int]], 
    domain2_residues: List[Tuple[str, int]], 
    engine: str = "cpptraj"
) -> np.ndarray:
    """Calculate the Domain-Domain Interaction Index (DSI) for a trajectory.

    DSI is a measure of distance between two domains, defined as:
    DSI = d(com1, com2) - (Rg1 + Rg2)
    where d(com1, com2) is the distance between the centers of mass of the
    two domains, and Rg1 and Rg2 are their respective radii of gyration.

    Args:
        ensemble: A StructureEnsemble object containing topology and trajectory.
        domain1_residues: A list of residue keys (chain_id, residue_idx)
                          for the first domain. If two keys for the same
                          chain are provided, they are treated as the
                          start and end of a continuous residue range.
        domain2_residues: A list of residue keys for the second domain.
        engine: The engine used for the DSI calculation. Default is "cpptraj".

    Returns:
        np.ndarray: DSI values for each frame in the trajectory.

    Raises:
        ValueError: If engine is not supported or if residue lists are invalid.
        TypeError: If ensemble is not a StructureEnsemble.

    Example:
        >>> from enzy_htp.structure import StructureEnsemble
        >>> from enzy_htp.analysis.dsi import dsi
        >>> 
        >>> # Define domain residues
        >>> domain1_res = [("A", 1), ("A", 50)]  # Chain A, residues 1-50
        >>> domain2_res = [("A", 100), ("A", 150)]  # Chain A, residues 100-150
        >>> 
        >>> # Calculate DSI
        >>> dsi_values = dsi(ensemble, domain1_res, domain2_res)
    """
    # Validate inputs
    if not isinstance(ensemble, StructureEnsemble):
        _LOGGER.error(f"ensemble must be a StructureEnsemble object. found: {type(ensemble)}")
        raise TypeError("ensemble must be a StructureEnsemble object")
    
    if engine not in DSI_METHODS:
        _LOGGER.error(f"engine ({engine}) not supported. Supported: {list(DSI_METHODS.keys())}")
        raise ValueError(f"Unsupported engine: {engine}")
    
    if not domain1_residues or not domain2_residues:
        _LOGGER.error("Both domain1_residues and domain2_residues must be non-empty")
        raise ValueError("Domain residue lists cannot be empty")
    
    # Expand residue ranges if needed
    domain1_expanded = _expand_residue_ranges(domain1_residues)
    domain2_expanded = _expand_residue_ranges(domain2_residues)
    
    # Dispatch to the appropriate method
    result = DSI_METHODS[engine](ensemble, domain1_expanded, domain2_expanded)
    
    return result


def _expand_residue_ranges(residues: List[Tuple[str, int]]) -> List[Tuple[str, int]]:
    """Expand residue ranges to individual residues when only 2 residues are provided.
    
    If two consecutive residues for the same chain are provided, they are treated
    as the start and end of a continuous residue range.
    
    Args:
        residues: List of (chain_id, residue_idx) tuples
        
    Returns:
        List of expanded individual residues
    """
    if len(residues) == 2:
        if residues[0][0] != residues[1][0]:
            # Different chains, treat as individual residues
            _LOGGER.warning(f"Two residues from different chains provided: {residues}. Treating as individual residues.")
            expanded = residues
        else:
            # Same chain, treat as a range
            start = residues[0][1]
            end = residues[1][1]
            expanded = [(residues[0][0], i) for i in range(start, end + 1)]
        return expanded
    else:
        return residues

# Method dispatch dictionary
DSI_METHODS: Dict[str, Callable] = {
    "cpptraj": interface.amber.calculate_dsi_metrics,
}
"""Submodule contains code for calculations the pKa of residues in the enzyme.
+ residue_pka()
    calculate the pKa of residues in a protein from a Structure or an ensemble
    of Structure()s

Author: QZ Shao <shaoqz@icloud.com>
Author: Robbie Ge <robbie.ge@vanderbilt.edu>

Date: 2025-08-07
"""
from typing import List, Dict, Union, Callable, Tuple, Optional

from enzy_htp import interface, config as eh_config
from enzy_htp.structure import Structure, StructureEnsemble, Residue
from enzy_htp.core import _LOGGER
from enzy_htp.preparation.clean import remove_solvent


def _parse_chain_res_code(res_str: str) -> Optional[Tuple[str, int]]:
    """Parse a residue code in 'chain.number' format into (chain_id, residue_number)."""
    parts = res_str.split('.')
    if len(parts) != 2:
        _LOGGER.warning(f"Invalid format for residue string '{res_str}'. Expected 'chain.number'")
        return None
    chain_id, num_str = parts
    try:
        return chain_id, int(num_str)
    except ValueError:
        _LOGGER.warning(f"Could not parse residue number from '{res_str}'")
        return None


def _parse_numeric_res_code(res_str: str) -> Optional[Tuple[str, int]]:
    """Parse a numeric residue code, assuming chain 'A'."""
    try:
        num = int(res_str)
        _LOGGER.warning(f"No chain specified for residue '{res_str}', assuming chain A")
        return "A", num
    except ValueError:
        _LOGGER.warning(f"Could not parse residue number from '{res_str}'")
        return None


def _parse_residue_key(res: Union[str, Residue]) -> Optional[Tuple[str, int]]:
    """Parse a residue identifier (str or Residue) into a (chain_id, residue_number) tuple."""
    if isinstance(res, str):
        # Try parsing 'chain.number' format
        parsed = _parse_chain_res_code(res)
        if parsed is not None:
            return parsed
        # Try numeric-only format
        numeric = _parse_numeric_res_code(res)
        if numeric is not None:
            return numeric
    elif isinstance(res, Residue):
        return res.key()
    else:
        _LOGGER.warning(f"Unsupported target_residues type: {type(res)}")
    return None


def _filter_target_residues(
        all_pka: Dict[Tuple[str, int], float], 
        target_residues: Union[List[str], List[Residue], None]
    ) -> Dict[Tuple[str, int], float]:
    """Helper function to filter pKa results for target residues using chain_id + res_num for robust alignment."""
    if target_residues is None:
        return all_pka
    
    # Convert target_residues to (chain_id, res_num) tuples
    target_keys = []
    for res in target_residues:
        key = _parse_residue_key(res)
        if key is not None:
            target_keys.append(key)
    
    # Filter results using (chain_id, res_num) keys for robust alignment
    return {key: pka for key, pka in all_pka.items() if key in target_keys}


def residue_pka(
        stru: Union[Structure, StructureEnsemble],
        target_residues: Union[List[Residue], List[str], None] = None,
        method: str = "propka",
        work_dir: str = None,
        remove_solvents: bool = True,
        **kwargs,
) -> Union[Dict[Tuple[str, int], float], List[Dict[Tuple[str, int], float]]]:
    """Calculate the pKa of target residues in a protein structure.
    Science API function for calculating pKa values using structure-based methods.
    The methods in this function are structure-based, meaning that a pKa value
    is calculated for a residue per structure. If the input is a StructureEnsemble,
    the pKa values will be calculated for each structure in the ensemble.
    On the contrary, there are other methods that calculate pKa values based on constant 
    pH MD simulations, which are not implemented in this function. 
    (there will be `residue_pka_md` in the future)

    Args:
        stru:
            The target protein as a Structure() or StructureEnsemble().
        target_residues:
            A list of residue keys (e.g.: "A.100") or Residue objects. If None,
            all ionizable residues will be analyzed.
        method:
            The algorithm for the pKa calculation. (see Details)
        work_dir:
            The working directory for the calculation. If None, uses the enzyhtp system.SCRATCH directory.
        remove_solvents:
            Whether to remove solvent molecules before calculation. Default True 
            since PROPKA doesn't benefit from solvents and they can cause formatting issues.

    Returns:
        A dictionary mapping (chain_id, residue_number) tuples to their pKa values, 
        or a list of such dictionaries for a StructureEnsemble.

    Details:
        Available implementations:
        - "propka": Uses PROPKA to calculate pKa values.
    """
    # Set default work_dir to SCRATCH if not provided
    if work_dir is None:
        work_dir = eh_config.system.SCRATCH_DIR
    
    if method not in PKA_METHODS:
        err_msg = f"Method '{method}' not supported. Supported methods: {list(PKA_METHODS.keys())}"
        _LOGGER.error(err_msg)
        raise ValueError(err_msg)

    if isinstance(stru, Structure):
        # Remove solvents if requested
        if remove_solvents:
            stru_clean = remove_solvent(stru, in_place=False)
        else:
            stru_clean = stru
        
        all_pka = PKA_METHODS[method](stru=stru_clean, work_dir=work_dir, **kwargs)
        return _filter_target_residues(all_pka, target_residues)

    elif isinstance(stru, StructureEnsemble):
        results = []
        for s in stru.structures(remove_solvent=remove_solvents):
            all_pka = PKA_METHODS[method](stru=s, work_dir=work_dir, **kwargs)
            results.append(_filter_target_residues(all_pka, target_residues))
        return results

    else:
        err_msg = f"stru can only be a Structure() or StructureEnsemble(). Found: {type(stru)}"
        _LOGGER.error(err_msg)
        raise TypeError(err_msg)

PKA_METHODS: Dict[str, Callable] = {
    "propka": interface.propka.get_residue_pka_from_stru,
}
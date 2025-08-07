"""Submodule contains code for calculations the pKa of residues in the enzyme.
+ residue_pka()
    calculate the pKa of residues in a protein from a Structure or an ensemble
    of Structure()s

Author: QZ Shao <shaoqz@icloud.com>
Author: Robbie Ge

Date: 2025-08-07
"""
from typing import List, Dict, Union, Callable, Tuple

from enzy_htp import interface, config as eh_config
from enzy_htp.structure import Structure, StructureEnsemble, Residue
from enzy_htp.core import _LOGGER
from enzy_htp.preparation.clean import remove_solvent


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
        if isinstance(res, str):
            # Format: "A.100" (chain.residue_number)
            if "." in res:
                parts = res.split(".")
                if len(parts) == 2:
                    chain_id = parts[0]
                    try:
                        res_num = int(parts[1])
                        target_keys.append((chain_id, res_num))
                    except ValueError:
                        _LOGGER.warning(f"Could not parse residue number from '{res}'")
                else:
                    _LOGGER.warning(f"Invalid format for residue string '{res}'. Expected 'chain.number'")
            else:
                # Just a number, assume chain A 
                try:
                    res_num = int(res)
                    target_keys.append(("A", res_num))
                    _LOGGER.warning(f"No chain specified for residue '{res}', assuming chain A")
                except ValueError:
                    _LOGGER.warning(f"Could not parse residue number from '{res}'")
        elif isinstance(res, Residue):
            target_keys.append(res.key())
        else:
            _LOGGER.warning(f"Unsupported target_residues type: {type(res)}")
    
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
            The working directory for the calculation. If None, uses the system SCRATCH directory.
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
        _LOGGER.error(f"Method '{method}' not supported. Supported methods: {list(PKA_METHODS.keys())}")
        raise ValueError

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
        _LOGGER.error(f"stru can only be a Structure() or StructureEnsemble(). Found: {type(stru)}")
        raise TypeError

PKA_METHODS: Dict[str, Callable] = {
    "propka": interface.propka.get_residue_pka_from_stru,
}
"""Submodule contains code for calculations the pKa of residues in the enzyme.
+ residue_pka()
    calculate the pKa of residues in a protein from a Structure or an ensemble
    of Structure()s

Author: QZ Shao <shaoqz@icloud.com>
Author: Robbie Ge

Date: 2025-08-07
"""
from typing import List, Dict, Union, Callable
from functools import partial
import os

from enzy_htp import interface
from enzy_htp.structure import Structure, StructureEnsemble, Residue
from enzy_htp.core import _LOGGER

def residue_pka(
        stru: Union[Structure, StructureEnsemble],
        target_residues: Union[List[Residue], List[str], None] = None,
        method: str = "propka",
        work_dir: str = "./pka",
        keep_in_file: bool = False,
        **kwargs,
) -> Union[Dict[int, float], List[Dict[int, float]]]:
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
            The working directory for the calculation.
        keep_in_file:
            Whether to keep the input files of the calculation.

    Returns:
        A dictionary of residue numbers and their pKa values, or a list of
        such dictionaries for a StructureEnsemble.

    Details:
        Available implementations:
        - "propka": Uses PROPKA to calculate pKa values.
    """
    if method not in PKA_METHODS:
        _LOGGER.error(f"Method '{method}' not supported. Supported methods: {list(PKA_METHODS.keys())}")
        raise ValueError

    def _filter_target_residues(all_pka: Dict[int, float], target_residues: Union[List[str], List[Residue], None]):
        """Helper function to filter pKa results for target residues."""
        if target_residues is None:
            return all_pka
        
        # Convert target_residues to residue numbers
        target_res_nums = []
        for res in target_residues:
            if isinstance(res, str):
                # Assume format like "A.100" or just "100" 
                if "." in res:
                    parts = res.split(".")
                    try:
                        target_res_nums.append(int(parts[-1]))
                    except ValueError:
                        _LOGGER.warning(f"Could not parse residue number from '{res}'")
                else:
                    try:
                        target_res_nums.append(int(res))
                    except ValueError:
                        _LOGGER.warning(f"Could not parse residue number from '{res}'")
            elif isinstance(res, Residue):
                target_res_nums.append(res.idx)
            else:
                _LOGGER.warning(f"Unsupported target_residues type: {type(res)}")
        
        # Filter results
        return {res_num: pka for res_num, pka in all_pka.items() if res_num in target_res_nums}

    if isinstance(stru, Structure):
        all_pka = PKA_METHODS[method](stru=stru, work_dir=work_dir, **kwargs)
        return _filter_target_residues(all_pka, target_residues)

    elif isinstance(stru, StructureEnsemble):
        results = []
        for s in stru.structures(remove_solvent=True):
            all_pka = PKA_METHODS[method](stru=s, work_dir=work_dir, **kwargs)
            results.append(_filter_target_residues(all_pka, target_residues))
        return results

    else:
        _LOGGER.error(f"stru can only be a Structure() or StructureEnsemble(). Found: {type(stru)}")
        raise TypeError

PKA_METHODS: Dict[str, Callable] = {
    "propka": interface.propka.get_residue_pka_from_stru,
}
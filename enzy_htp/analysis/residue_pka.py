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
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.structure.structure_io import PDBParser

def residue_pka(
        stru: Union[Structure, StructureEnsemble],
        target_residues: Union[List[Residue], List[str]],
        method: str = "propka",
        work_dir: str = "./pka",
        keep_in_file: bool = False,
        **kwargs,
) -> Union[Dict[int, float], List[Dict[int, float]]]:
    """Calculate the pKa of target residues in a protein structure.
    The methods in this function are structure-based, meaning that the a pKa value
    is calculated for a residue per structure. If the input is a StructureEnsemble,
    the pKa values will be calculated for each structure in the ensemble.
    On the contrary, there are other methods that calculate pKa values based on constant 
    pH MD simulations, which are not implemented in this function. 
    (there will be `residue_pka_md` in the future)

    Args:
        stru:
            The target protein as a Structure() or StructureEnsemble().
        target_residues:
            A list of residue keys (e.g.: "A.100") or Residue objects.
        method:
            The algorithm for the pKa calculation. (see Details)
        cluster_job_config:
            The config for ClusterJob. (enables running via a ClusterJob/ARMer)
        job_check_period:
            The time cycle for updating job state change (Unit: s).
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

    # Parse target residues
    # TODO use the get method
    target_res_list = None # place holder

    if isinstance(stru, Structure):
        all_pka = PKA_METHODS[method](stru=stru, work_dir=work_dir, **kwargs)
        return {res: all_pka.get(res) for res in target_res_list}

    elif isinstance(stru, StructureEnsemble):
        results = []
        for s in stru.structures(remove_solvent=True):
            all_pka = PKA_METHODS[method](stru=s, work_dir=work_dir, **kwargs)
            results.append({res: all_pka.get(res) for res in target_res_list})
        return results

    else:
        _LOGGER.error(f"stru can only be a Structure() or StructureEnsemble(). Found: {type(stru)}")
        raise TypeError

PKA_METHODS: Dict[str, Callable] = {
    "propka": interface.propka.get_residue_pka_from_stru,
}
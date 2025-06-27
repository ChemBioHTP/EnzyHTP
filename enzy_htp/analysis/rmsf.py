"""Submodule contains code for calculations the RMSF values from a structure ensemble
+ rmsf()
    Calculate the RMSF value of a StructureEnsemble instance.

Author: Qianzhen Shao <shaoqz@icloud.com>

Date: 2025-4-24
"""
from typing import Dict

from enzy_htp.preparation import remove_solvent
from enzy_htp import interface as eh_interface
from enzy_htp.structure import StruSelection
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru

def rmsf(stru_esm: StructureEnsemble, region_pattern: str = "polymer and (not elem H)",
        by_residue: bool = True, ignore_solvent: bool = True) -> Dict[str, float]:
    """Calculate the RMSF values of each atoms in the region_pattern of a StructureEnsemble
    instance. RMSF is calculated based on the average structure.
    Powered by the atomicfluct from Cpptraj from AmberTools for now.

    Args:
        stru_esm: 
            A conformational ensemble of a structure.
        region_pattern: 
            A pymol-style selection pattern that defines atoms of the RMSF calculation.
        by_residue:
            control if return values are grouped by residues or not. If True, calculate the RMSF by residues
            by taking the square-mean-root of the RMSF of each composing atoms.
        ignore_solvent:
            control of solvent are removed before the calculation. True by default to speed up the calculation
            if solvent is the interest of study, set this to False.

    Returns:
        A dictionary that map a EnzyHTP get pattern (<chain_name>.<residue_index>.<atom_name>) to the RMSF value.
        Example: {"A.1.CA" : 0.1} or {"A.2" : 1.1}
    """
    if ignore_solvent:
        stru = remove_solvent(stru_esm.structure_0)
    else:
        stru = stru_esm.structure_0
    stru_sele: StruSelection = select_stru(stru, pattern=region_pattern)
    return eh_interface.amber.get_rmsf(stru_esm=stru_esm, stru_selection=stru_sele, by_residue=by_residue)

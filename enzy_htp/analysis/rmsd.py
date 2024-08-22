#! python3
# -*- encoding: utf-8 -*-
'''
@File    :   rmsd.py
@Created :   2024/07/12 18:45
@Author  :   Zhong, Yinjie
@Version :   2.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
from typing import List

# Here put EnzyHTP imports.
from enzy_htp import interface, _LOGGER
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs
from enzy_htp.structure import Residue
from enzy_htp.structure.structure_ensemble import StructureEnsemble

def compose_mask_pattern(mask_region: List[Residue] = list(), ca_only: bool = True) -> str:
    """Compose the mask pattern from the given mask region.
    
    Args:
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    mask_pattern = str()
    if (mask_region):
        mask_pattern = f"resi {'+'.join([str(residue.idx) for residue in mask_region])}"
    else:
        mask_pattern = "all"

    if (ca_only):
        mask_pattern += " and name CA"
    return mask_pattern

def rmsd_of_structure(structure_ensemble: StructureEnsemble, include_ligand: bool = True, ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result of the whole structure.
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    mask_region = structure_ensemble.structure_0.residues
    if (not include_ligand):
        mask_region = list(filter(lambda residue: (not residue.is_noncanonical()), mask_region))
    return rmsd_of_region(structure_ensemble, mask_region, ca_only)

def rmsd_of_region(structure_ensemble: StructureEnsemble, mask_region: List[Residue] = list(), ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result with a mask region (a list of residues).
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    mask_pattern = compose_mask_pattern(mask_region=mask_region, ca_only=ca_only)
    return rmsd_with_pattern(structure_ensemble, mask_pattern)

def rmsd_with_pattern(structure_ensemble: StructureEnsemble, mask_pattern: str) -> float:
    """
    mvp function for RMSD calculation
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        mask_pattern (str): A pymol-formatted selection string which defines the region for calculating RMSD value.
    """
    return interface.pymol.get_rmsd(structure_ensemble=structure_ensemble, mask_pattern=mask_pattern)

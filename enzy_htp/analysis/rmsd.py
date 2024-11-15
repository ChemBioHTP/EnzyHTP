"""Submodule contains code for calculations the RMSD value of a StructureEnsemble instance.
The calculation can be performed with Amber interface.
+ rmsd_calculation()
    Calculate the RMSD value of a StructureEnsemble instance with specified mask pattern.
+ rmsd_of_region()
    Get RMSD value from MD simulation result with a mask region (a list of residues).
+ rmsd_of_structure()
    Get RMSD value from MD simulation result of the whole structure.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>

Date: 2024-11-15
"""
# Here put the import lib.
from typing import List, Dict, Union

# Here put EnzyHTP imports.
from enzy_htp import interface as eh_interface
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs
from enzy_htp.preparation import remove_solvent
from enzy_htp.structure import Residue, Structure
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru

def compose_pymol_pattern(mask_region: List[Residue] = list(), ca_only: bool = True) -> str:
    """Compose the given list of Residue instances to PyMOL-format pattern.
    
    Args:
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms are included. Default True.
    
    Returns:
        region_pattern (str): A pymol-formatted selection string which defines the region for calculating RMSD value.
    """
    region_pattern = str()

    if (mask_region):
        region_pattern = f"resi {'+'.join([str(residue.idx) for residue in mask_region])}"
    else:
        region_pattern = "all"
    if (ca_only):
        region_pattern += " and name CA"
    else:
        region_pattern += " and (not elem H)"
    return region_pattern

def rmsd_calculation(stru_esm: StructureEnsemble, region_pattern: str = "all") -> float:
    """Calculate the RMSD value of a StructureEnsemble instance with specified mask pattern.
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        region_pattern (str): A pymol-formatted selection string which defines the region for calculating RMSD value.
    """
    stru = remove_solvent(stru_esm.structure_0, in_place=False)
    stru_sele = select_stru(stru, pattern=region_pattern)
    return eh_interface.amber.get_rmsd(traj_path=stru_esm.coordinate_list, prmtop_path=stru_esm.topology_source_file, mask_selection=stru_sele)

def rmsd_of_region(stru_esm: StructureEnsemble, mask_region: List[Residue] = list(), ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result with a mask region (a list of residues).
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    region_pattern = compose_pymol_pattern(mask_region=mask_region, ca_only=ca_only)
    return rmsd_calculation(
        stru_esm=stru_esm, 
        region_pattern=region_pattern
    )

def rmsd_of_structure(stru_esm: StructureEnsemble, include_ligand: bool = True, ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result of the whole structure. Solvents are not included.
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    stru = remove_solvent(stru_esm.structure_0, in_place=False)
    mask_region = stru.residues
    if (not include_ligand):
        mask_region = list(filter(lambda residue: (not residue.is_noncanonical()), mask_region))
    return rmsd_of_region(
        stru_esm=stru_esm, 
        mask_region=mask_region, 
        ca_only=ca_only
    )

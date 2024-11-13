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
from typing import List, Dict, Union

# Here put EnzyHTP imports.
from enzy_htp import interface as eh_interface, _LOGGER
from enzy_htp import config as eh_config
from enzy_htp._interface.pymol_interface import PyMolInterface
from enzy_htp._interface.amber_interface import AmberInterface
from enzy_htp.core import file_system as fs
from enzy_htp.structure import Residue, Structure
from enzy_htp.structure.structure_ensemble import StructureEnsemble

ACCEPTED_INTERFACES: Dict[str, Union[PyMolInterface, AmberInterface]] = {
    "pymol": eh_interface.pymol,
    "amber": eh_interface.amber,
}

def compose_mask_pattern(mask_region: List[Residue] = list(), interface: str = "amber", ca_only: bool = True) -> str:
    """Compose the mask pattern from the given mask region for Amber/CppTraj or PyMOL interface.
    
    Args:
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        interface (str, optional): The interface for calculation. Default amber.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    interface = interface.lower()
    mask_pattern = str()
    
    if (interface == "amber"):
        if (mask_region):
            mask_pattern = f":{','.join([str(residue.idx) for residue in mask_region])}"
        else:
            mask_pattern = ""
        if (ca_only):
            if mask_pattern:
                mask_pattern += "&"
            mask_pattern += "@CA"
            
    elif (interface == "pymol"):
        if (mask_region):
            mask_pattern = f"resi {'+'.join([str(residue.idx) for residue in mask_region])}"
        else:
            mask_pattern = "all"
        if (ca_only):
            mask_pattern += " and name CA"

    elif (interface not in ACCEPTED_INTERFACES.keys()):
        raise ValueError(f"Unacceptable interface value. Please choose from {' and '.join(ACCEPTED_INTERFACES.keys())}")
    return mask_pattern

def rmsd_of_structure(structure_ensemble: StructureEnsemble, interface: str = "amber", reference_structure: Structure = None, include_ligand: bool = True, ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result of the whole structure.
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        interface (str, optional): The interface for calculation. Default amber.
        reference_structure (Structure, optional): The reference structure, which should have the same sequences and ligangs as the structure_ensemble. Default: Average structure.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    interface = interface.lower()
    mask_region = structure_ensemble.structure_0.residues
    if (not include_ligand):
        mask_region = list(filter(lambda residue: (not residue.is_noncanonical()), mask_region))
    return rmsd_of_region(
        structure_ensemble=structure_ensemble, 
        interface=interface,
        reference_structure=reference_structure, 
        mask_region=mask_region, 
        ca_only=ca_only
    )

def rmsd_of_region(structure_ensemble: StructureEnsemble, interface: str = "amber", reference_structure: Structure = None, mask_region: List[Residue] = list(), ca_only: bool = True) -> float:
    """Get RMSD value from MD simulation result with a mask region (a list of residues).
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        interface (str, optional): The interface for calculation. Default amber.
        reference_structure (Structure, optional): The reference structure, which should have the same sequences and ligangs as the structure_ensemble. Default: Average structure.
        mask_region (List[Residue], optional): A list of Residues to calculation RMSD value. If the list is empty, then all residues are selected.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    interface = interface.lower()
    mask_pattern = compose_mask_pattern(mask_region=mask_region, interface=interface, ca_only=ca_only)
    return rmsd_with_pattern(
        structure_ensemble=structure_ensemble, 
        interface=interface,
        reference_structure=reference_structure, 
        mask_pattern=mask_pattern
    )

def rmsd_with_pattern(structure_ensemble: StructureEnsemble, interface: str = "amber", reference_structure: Structure = None, mask_pattern: str = str()) -> float:
    """
    mvp function for RMSD calculation
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        interface (str, optional): The interface for calculation. Default amber.
        reference_structure (Structure, optional): The reference structure, which should have the same sequences and ligangs as the structure_ensemble. Default: Average structure.
        mask_pattern (str): A pymol-formatted selection string which defines the region for calculating RMSD value.
    """
    interface = interface.lower()
    if (interface == "amber"):
        return ACCEPTED_INTERFACES[interface].get_rmsd(traj_path=structure_ensemble.coordinate_list, prmtop_path=structure_ensemble.topology_source_file, mask_pattern=mask_pattern)
    elif (interface == "pymol"):
        return ACCEPTED_INTERFACES[interface].get_rmsd(structure_ensemble=structure_ensemble, reference_structure=reference_structure, mask_pattern=mask_pattern)

    elif (interface not in ACCEPTED_INTERFACES.keys()):
        raise KeyError(f"Unacceptable interface value. Please choose from {' and '.join(ACCEPTED_INTERFACES.keys())}")

#! python3
# -*- encoding: utf-8 -*-
'''
@File    :   rmsd.py
@Created :   2024/07/12 18:45
@Author  :   Zhong, Yinjie
@Version :   1.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
import os
from subprocess import CalledProcessError, run
from pandas import read_csv
from typing import List

# Here put EnzyHTP imports.
from enzy_htp import interface, _LOGGER
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs
from enzy_htp.structure import Residue
from enzy_htp.structure.structure_ensemble import StructureEnsemble

def compose_mask_pattern(mask_region: List[Residue], ca_only: bool = False) -> str:
    """Compose the mask pattern from the given mask region.
    
    Args:
        mask_region (List[Residue]): A list of Residues to calculation RMSD value.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default False.
    """
    mask_pattern = str()
    if (ca_only):
        mask_pattern = f":{','.join([str(residue.idx) for residue in mask_region])}@CA"
    else:
        mask_pattern = f":{','.join([str(residue.idx) for residue in mask_region])}&!@H="
    return mask_pattern

def rmsd_of_structure(structure_ensemble: StructureEnsemble, include_ligand: bool = True, ca_only: bool = False) -> float:
    """Get RMSD value from MD simulation result of the whole structure.
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default False.
    """
    mask_region = structure_ensemble.structure_0.residues
    if (not include_ligand):
        mask_region = list(filter(lambda residue: (not residue.is_noncanonical()), mask_region))
    return rmsd_of_region(structure_ensemble, mask_region, ca_only)

def rmsd_of_region(structure_ensemble: StructureEnsemble, mask_region: List[Residue], ca_only: bool = False) -> float:
    """Get RMSD value from MD simulation result with a mask region (a list of residues).
    
    Args:
        structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        mask_region (List[Residue]): A list of Residues to calculation RMSD value.
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
    tmp_rmsd_in = f"{eh_config['system.SCRATCH_DIR']}/cpptraj_rmsd.in"
    tmp_rmsd_out = f"{eh_config['system.SCRATCH_DIR']}/rmsd.agr"

    prmtop_path = structure_ensemble._topology
    traj_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/cpptraj.in")
    traj_path = structure_ensemble.coordinate_list

    with open(tmp_rmsd_in, "w") as fobj:
        fobj.write(f"""parm {prmtop_path}
trajin {traj_path}
autoimage
rmsd {mask_pattern} first mass
average crdset AVE {mask_pattern}
run
autoimage
rmsd {mask_pattern} ref AVE * out {tmp_rmsd_out} mass
run
""")
    try:
        run(f"cpptraj -i {tmp_rmsd_in}",
        check=True, text=True, shell=True, capture_output=True)
    except CalledProcessError as e:
        _LOGGER.error(f"{e.stdout}{os.linesep}{e.stderr}")
        raise e
    result = get_cpptraj_rmsd_result(tmp_rmsd_out)
    fs.safe_rm(tmp_rmsd_in)
    fs.safe_rm(tmp_rmsd_out)
    return result

def get_cpptraj_rmsd_result(result_path: str) -> float:
    """
    mvp function for RMSD calculation
    """
    result_df = read_csv(result_path, delim_whitespace=True)
    return result_df.iloc[:, 1].mean()
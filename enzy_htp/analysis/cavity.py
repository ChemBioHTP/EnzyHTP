"""Submodule contains code for cavity calculation of a Structure or StructureEnsemble instance.
The calculation is performed with Mole2 interface.
+ identify_stru_cavities()
    Identify the cavity of a Structure instance with specified configs.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>

Date: 2024-11-15
"""
from __future__ import annotations
from os import path
from typing import Callable, Dict, List, Tuple

import numpy as np
import numpy.typing as npt
import pyvista as pv

from enzy_htp import interface, _LOGGER, PDBParser
from enzy_htp.structure import Structure, StructureEnsemble, Residue, Cavity
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs

sp = PDBParser()
CAVITY_IDENTIFICATION_METHODS: Dict[str, Callable[..., List[Cavity]]] = {
    "mole2": interface.mole2.identify_cavities
}

def identify_stru_cavities(stru: Structure,
        work_dir: str = None,
        engine: str = "mole2",
        **kwargs
    ) -> List[Cavity]:
    """Identifies cavities in a Structure instance. 
    Currently available engine:
    * mole2;

    Args:
        stru (Structure): The structure instance to detect cavities from.
        work_dir (str, optional): Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
        engine (str, optional): The engine to use for cavity identification. Defaults to "mole2" (The only available one at present).
        **kwargs: Engine-specific parameters in keyword arguments. 

    Returns:
        cavities (List[Cavity]): A list of Cavity objects.

    Details:
        * mole2 specific arguments:
            - non_active_residues (List[Residue], optional): Residues that should be skipped.
            - probe (float, optional): Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
            - inner (float, optional): Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
            - mesh_density (float, optional): Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
            - ignore_hetatm (bool, optional): TODO (CJ)
            - use_mono (bool, optional): Indicate if mono need to be used during run time. Defaults to true.

    """
    if work_dir is None:
        work_dir = eh_config['system.SCRATCH_DIR']

    cavities = CAVITY_IDENTIFICATION_METHODS[engine](stru=stru, work_dir=work_dir, **kwargs)
    
    return cavities

def _choose_cavity(cavity_list: List[Cavity], 
        composing_residues: List[Residue] = None, 
        contain_ligand: str = None, 
        target_cavity: Cavity = None
    ) -> Tuple[Cavity, int]:
    """
    Selects the most similar cavity (to the focus cavity) from the `cavity_list` based on residue key overlap.
    The cavity with most residue key overlap is selected.
    
    Args:
        cavity_list: List of Cavity objects to compare. All cavities should be in the same Structure instance.
        composing_residues (List[Residue], optional): A list of Residue instances composing the cavity. Mutually exclusive with `contain_ligand` or `target_cavity`.
        contain_ligand (str, optional): The PyMOL-formatted pattern selecting the ligand in the cavity. Mutually exclusive with `target_cavity` or `composing_residues`.
        target_cavity (Cavity, optional): The target cavity to track throughout the ensemble (by most overlapped residues). Mutually exclusive with `composing_residues` or `contain_ligand`.
    
    Returns:
        result (Tuple[Cavity, int]): 
            * Cavity: The Cavity object from cavity_list with highest residue key overlap with focus_cavity
            * int: max_overlap value.
    """
    if sum(x is not None for x in (composing_residues, contain_ligand, target_cavity)) != 1:
        _LOGGER.error("The `composing_residues`, `contain_ligand` and `target_cavity` are mutually exclusive to each other.")
        raise ValueError()
    
    if (contain_ligand):
        stru = cavity_list[0].stru
        ligand_selection = select_stru(stru=stru, pattern=contain_ligand)
        ligand_atom_points = [np.array(atom.coord) for atom in ligand_selection.atoms]
        if (len(ligand_atom_points) == 0):  # If nothing selected.
            return None, 0
        cavity_point_dict = dict()  # A dict recording how many atoms are contained by each cavity.
        for cavity in cavity_list:
            contain_point_list = cavity.contains_array(ligand_atom_points)
            contain_point_count = sum(1 for point in contain_point_list if point == True)   # Number of atoms in the cavity.
            if (contain_point_count > 0):
                cavity_point_dict[cavity] = contain_point_count
            continue
        _LOGGER.info(cavity_point_dict)
        if (len(cavity_point_dict.keys()) > 0):
            # Return the cavity containing most atoms of the ligand.
            selected_cavity = cavity_list[max(cavity_point_dict, key=cavity_point_dict.get)]
            return selected_cavity, 1
        else:
            # If none of the cavities containing any atom of the ligand, return default value.
            return None, 0
    else:
        focus_residue_keys = set()
        cavity_similarity = dict()
        if (target_cavity):
            focus_residue_keys = set(resi.key() for resi in (target_cavity.boundary_residues + target_cavity.inner_residues))
        if (composing_residues):
            focus_residue_keys = set(resi.key() for resi in composing_residues)
        for i, cavity in enumerate(cavity_list):
            cavity_residue_keys = set(resi.key() for resi in (cavity.boundary_residues + cavity.inner_residues))
            cavity_similarity[i] = len(cavity_residue_keys.intersection(focus_residue_keys))
        return cavity_list[max(cavity_similarity, key=cavity_similarity.get)], max(cavity_similarity.values())

def ensemble_cavity_volumes(
        stru_esm: StructureEnsemble, 
        composing_residues: List[Residue] = None, 
        contain_ligand: str = None, 
        target_cavity: Cavity = None, 
        frame_0_based: bool = True,
        work_dir: str = None,
        engine: str = "mole2",
        **kwargs) -> List[float]:
    """Caculate cavity volumes from a StructureEnsemble instance, which is done by each frame.
    The cavity should be selected by a list of Residue instances, a ligand string (PyMOL formatted), or a Cavity instance.
    
    Args:
        stru_esm (StructureEnsemble): The StructureEnsemble instance to calculate cavities from.
        composing_residues (List[Residue], optional): A list of Residue instances composing the cavity. Mutually exclusive with `contain_ligand` or `target_cavity`.
        contain_ligand (str, optional): The PyMOL-formatted pattern selecting the ligand in the cavity. Mutually exclusive with `target_cavity` or `composing_residues`.
        target_cavity (Cavity, optional): The target cavity to track throughout the ensemble (by most overlapped residues). Mutually exclusive with `composing_residues` or `contain_ligand`.
        frame_0_based (bool, optional): Indicate if the cavity selection input is formulated with the frame 0 structure to get target cavity.
        work_dir (str, optional): Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
        engine (str, optional): The engine to use for cavity identification.

    Returns:
        volumes (List[float]): The list of cavity volume value of each frame from the ensemble.
    """
    if sum(x is not None for x in (composing_residues, contain_ligand, target_cavity)) != 1:
        _LOGGER.error("The `composing_residues`, `contain_ligand` and `target_cavity` are mutually exclusive to each other.")
        raise ValueError()
    esm_cavities: List[Cavity] = list()
    structure_0 = stru_esm.structure_0

    confirmed_target_cavity = None
    if frame_0_based:   # Confirm the target cavity if `frame_0_based=True`.
        non_active_residues = []
        if (kwargs.get(contain_ligand)):
            ligand_selection = select_stru(stru=structure_0, pattern=kwargs.get(contain_ligand))
            non_active_residues = ligand_selection.involved_residues
        frame_0_cavities = identify_stru_cavities(stru=structure_0, 
            work_dir=work_dir, engine=engine, 
            non_active_residues=non_active_residues, **kwargs)
        confirmed_target_cavity, _ = _choose_cavity(cavity_list=frame_0_cavities, 
            composing_residues=composing_residues, contain_ligand=contain_ligand, target_cavity=target_cavity)
        if (confirmed_target_cavity is None):
            _LOGGER.error("Unable to identify target cavity from frame 0 structure.")
            raise ValueError()
        else:
            # If the target cavity is confirmed, we will track the target cavity instead of `composing_residues` or `contain_ligand`.
            composing_residues = None
            contain_ligand = None
    else:
        confirmed_target_cavity = target_cavity
    
    for stru_frame, _, _ in stru_esm.structures(remove_solvent=True):     # Iterate over the ensemble.
        non_active_residues = []
        if (kwargs.get(contain_ligand)):
            ligand_selection = select_stru(stru=structure_0, pattern=kwargs.get(contain_ligand))
            non_active_residues = ligand_selection.involved_residues
        frame_cavities = identify_stru_cavities(stru=stru_frame, work_dir=work_dir, engine=engine, 
            non_active_residues=non_active_residues, **kwargs)
        cavity, max_overlap = _choose_cavity(cavity_list=frame_cavities, 
            composing_residues=composing_residues, contain_ligand=contain_ligand, target_cavity=confirmed_target_cavity)
        esm_cavities.append(cavity)
        continue

    volumes = list()
    for cavity in esm_cavities:
        if (cavity is None):
            volumes.append(0.0)
        else:
            volumes.append(cavity.volume)
    return volumes

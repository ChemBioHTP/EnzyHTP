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
from enzy_htp import config as eh_config
from enzy_htp.structure.structure_selection import select_stru
import enzy_htp.structure.structure_operation as so
from enzy_htp.structure import Structure, StructureEnsemble, Residue, Cavity
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
    Selects the most relevant Cavity from a list according to exactly one of three modes:
        1. contain_ligand: find the cavity enclosing the most ligand atoms matching a PyMOL pattern;
        2. composing_residues: find the cavity sharing the most residue keys with a given residue list;
        3. target_cavity: find the cavity sharing the most residue keys with a reference cavity.

    Exactly one of `composing_residues`, `contain_ligand`, or `target_cavity` must be provided;
    otherwise a ValueError is raised.
    Args:
        cavity_list (List[Cavity]):
                List of candidate cavities (all must belong to the same Structure).
        composing_residues (List[Residue], optional):
                Residues defining the cavity of interest. Mutually exclusive with
                `contain_ligand` and `target_cavity`.
        contain_ligand (str, optional):
                PyMOL selection pattern for ligand atoms. Mutually exclusive with
                `composing_residues` and `target_cavity`.
        target_cavity (Cavity, optional):
                Reference cavity whose residues define the target site. Mutually exclusive with
                `composing_residues` and `contain_ligand`.
    Returns:
        Tuple[Cavity, int]:
            * The selected Cavity object from `cavity_list` with the highest residue-key overlap.
            * The confidence score from 0 to 1. (
                - for composing_residue: it is the ratio of shared residues to total residues in the cavity.
                - for contain_ligand: it is the ratio of ligand atoms contained in the cavity to total ligand atoms.
                - for target_cavity: it is a boolean indicating whether the cavity matches the target.
            )
    Raises:
        ValueError:
                If none or more than one of `composing_residues`, `contain_ligand`, or `target_cavity`
                is specified.

        Tuple[Cavity or None, int]:
                - In ligand-containment mode: the cavity enclosing the most ligand atoms (or None if
                    no atoms are found) and a flag 1 (if any cavity contains ligand) or 0.
                - In residue-overlap modes: the cavity with the highest count of overlapping residue
                    keys and the integer count of those shared residues.
    """
    if sum(x is not None for x in (composing_residues, contain_ligand, target_cavity)) != 1:
        err_msg = "The `composing_residues`, `contain_ligand` and `target_cavity` are mutually exclusive and at least one must be specified."
        _LOGGER.error(err_msg)
        raise ValueError(err_msg)

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
        _LOGGER.info(f"ligand containing cavities (cavity:number of atoms contained):{cavity_point_dict}")
        if (len(cavity_point_dict.keys()) > 0):
            # Return the cavity containing most atoms of the ligand.
            selected_cavity = max(cavity_point_dict, key=cavity_point_dict.get)
            n_atoms_in_cavity = cavity_point_dict[selected_cavity]
            ratio = n_atoms_in_cavity / len(ligand_atom_points)
            return selected_cavity, ratio
        else:
            # If none of the cavities containing any atom of the ligand, return default value.
            return None, 0
    else:
        focus_residue_keys = set()
        cavity_similarity = []
        if (target_cavity):
            focus_residue_keys = set(resi.key() for resi in (target_cavity.boundary_residues + target_cavity.inner_residues))
        if (composing_residues):
            focus_residue_keys = set(resi.key() for resi in composing_residues)
        for cavity in cavity_list:
            cavity_residue_keys = set(resi.key() for resi in (cavity.boundary_residues + cavity.inner_residues))
            cavity_similarity.append(len(cavity_residue_keys.intersection(focus_residue_keys)))
        return cavity_list[np.argmax(cavity_similarity)], max(cavity_similarity) / len(focus_residue_keys)

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
        err_msg = "The `composing_residues`, `contain_ligand` and `target_cavity` are mutually exclusive and at least one must be specified."
        _LOGGER.error(err_msg)
        raise ValueError(err_msg)
    esm_cavities: List[Cavity] = list()
    structure_0 = so.remove_solvent(stru_esm.structure_0) # dont need copy as this is lazy-generated property

    confirmed_target_cavity = None
    if frame_0_based:   # Confirm the target cavity if `frame_0_based=True`.
        non_active_residues = []
        if contain_ligand:
            ligand_selection = select_stru(stru=structure_0, pattern=contain_ligand)
            non_active_residues = ligand_selection.involved_residues
        frame_0_cavities = identify_stru_cavities(stru=structure_0, 
            work_dir=work_dir, engine=engine, 
            non_active_residues=non_active_residues, **kwargs)
        confirmed_target_cavity, confidence = _choose_cavity(cavity_list=frame_0_cavities, 
            composing_residues=composing_residues, contain_ligand=contain_ligand, target_cavity=target_cavity)
        if (confirmed_target_cavity is None):
            err_msg = "Unable to identify target cavity from frame 0 structure."
            _LOGGER.error(err_msg)
            raise ValueError(err_msg)
        else:
            # If the target cavity is confirmed, we will track the target cavity instead of `composing_residues` or `contain_ligand`.
            composing_residues = None
            contain_ligand = None
    else:
        confirmed_target_cavity = target_cavity
    
    for stru_frame, _, _ in stru_esm.structures(remove_solvent=True):     # Iterate over the ensemble.
        non_active_residues = []
        if contain_ligand:
            ligand_selection = select_stru(stru=structure_0, pattern=contain_ligand)
            non_active_residues = ligand_selection.involved_residues
        frame_cavities = identify_stru_cavities(stru=stru_frame, work_dir=work_dir, engine=engine, 
            non_active_residues=non_active_residues, **kwargs)
        cavity, confidence = _choose_cavity(cavity_list=frame_cavities, 
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

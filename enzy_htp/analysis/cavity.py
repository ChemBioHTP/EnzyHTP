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

from enzy_htp import interface, _LOGGER, PDBParser, Residue
from enzy_htp.structure import Structure, Residue, Cavity
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
        engine (str): The engine to use for cavity identification. Defaults to "mole2" (The only available one at present).
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

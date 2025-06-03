"""Submodule contains code for cavity calculation of a Structure or StructureEnsemble instance.
The calculation is performed with Mole2 interface.
+ identify_stru_cavities()
    Identify the cavity of a Structure instance with specified configs.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>

Date: 2024-11-15
"""
from __future__ import annotations
from os import path
from typing import List, Tuple

import numpy as np
import numpy.typing as npt
import pyvista as pv

from enzy_htp import interface, _LOGGER, PDBParser, Structure, Residue
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs
from enzy_htp._interface import Mole2Cavity
from enzy_htp.structure import StructureCavity

sp = PDBParser()
mole2 = interface.mole2

def identify_stru_cavities(stru: Structure,
        non_active_residues: List[Residue] = list(), 
        probe: float = None, 
        inner: float = None, 
        mesh_density: float = None,
        ignore_hetatm: bool = None,
        work_dir: str = None,
        use_mono: bool = True
    ) -> List[StructureCavity]:
    """Identifies cavities in a protein structure using the Mole2 software package. Client method that should be 
    called by users. Results are represented via Mole2Cavity objects that support basic geometry operations.

    Args:
        stru (Structure): The structure instance containing the cavity.
        non_active_residues (List[Residue], optional): Residues that should be skipped.
        probe (float, optional): Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
        inner (float, optional): Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
        mesh_density (float, optional): Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
        ignore_hetatm (bool, optional): TODO(CJ)
        work_dir (str, optional): Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
        use_mono (bool, optional): Does mono need to be used during run time? Defaults to true.

    Returns:
        A list() of Mole2Cavity ebjects.
    """

    if work_dir is None:
        work_dir = eh_config['system.SCRATCH_DIR']
    pdb_filepath = path.join(work_dir, "stru_cavity_temp.pdb")
    non_active_parts = [resi.key() for resi in non_active_residues]

    mole2cavities = mole2.identify_cavities(pdb_path=pdb_filepath, non_active_parts=non_active_parts,
        probe=probe, inner=inner, mesh_density=mesh_density, ignore_hetatm=ignore_hetatm,
        work_dir=work_dir, use_mono=use_mono)
    
    stru_cavities = StructureCavity.from_mole2_cavities(stru=stru, mole2cavities=mole2cavities)
    return stru_cavities

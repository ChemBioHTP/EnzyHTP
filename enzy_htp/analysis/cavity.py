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

from enzy_htp import interface, _LOGGER, PDBParser, Residue
from enzy_htp.structure import Structure
from enzy_htp import config as eh_config
from enzy_htp.core import file_system as fs
from enzy_htp._interface import Mole2Cavity

sp = PDBParser()
mole2 = interface.mole2

class Cavity(Mole2Cavity):
    """Class for the cavity of a Structure or StructureEnsemble instance.
    The calculation is performed with MOLE2 interface.
    """
    def __init__(self, stru: Structure, 
            points :List[npt.NDArray], mesh : pv.PolyData, probe:float, 
            inner: float, mesh_density:float, com: npt.NDArray, mole2_volume: float,            
            boundary_residue_keys: List[Tuple[str, int]], inner_residue_keys: List[Tuple[str, int]]):
        """Initialize a Cavity instance.
        
        Args:
            stru (Structure): The structure instance containing the cavity.
            points (List[npt.NDArray]): A list() of 3D vertices describing the cavity.
            mesh (pv.PolyData): The mesh representing the cavity established with pyvista.
            probe (float): Probe radius used during collection in Angstroms.
            inner (float): Inner radius used during collection in Angstroms.
            com (npt.NDArray): The center-of-mass of the mesh as a numpy array with format (x, y, z).
            mole2_volume (float): The volume of the cavity in A^3 calculated by Mole2 engine.
            boundary_residue_keys (List[Tuple[str, int]]): A list() of residue keys forming the cavity that are on the boundary of the structure.
            inner_residue_keys (List[Tuple[str, int]]): A list() of residue keys forming the cavity that are inside the structure.
        """
        self.stru = stru
        super().__init__(points=points, mesh=mesh, probe=probe, inner=inner, 
                mesh_density=mesh_density, com=com, mole2_volume=mole2_volume, 
                boundary_residue_keys=boundary_residue_keys, inner_residue_keys=inner_residue_keys)
        self.boundary_residues = [stru.find_residue_with_key(key) for key in self.boundary_residue_keys] if self.boundary_residue_keys else []
        self.inner_residues = [stru.find_residue_with_key(key) for key in self.inner_residue_keys] if self.inner_residue_keys else []
        return
    
    @classmethod
    def from_mole2_cavity(cls, stru: Structure, mole2cavity: Mole2Cavity) -> Cavity:
        """Create a Cavity instance from a Mole2Cavity instance.
        
        Args:
            stru (Structure): The structure instance containing the cavity.
            mole2cavity (Mole2Cavity): The Mole2Cavity Instance to initialize the Cavity instance.
        """
        cavity = cls(stru, points=mole2cavity.points(), 
                mesh=mole2cavity.mesh_, probe=mole2cavity.probe(), inner=mole2cavity.inner(), 
                mesh_density=mole2cavity.mesh_density(), com=mole2cavity.com_, mole2_volume=mole2cavity.mole2_volume, 
                boundary_residue_keys=mole2cavity.boundary_residue_keys, inner_residue_keys=mole2cavity.inner_residue_keys)
        return cavity
    
    @classmethod
    def from_mole2_cavities(cls, stru: Structure, mole2cavities: List[Mole2Cavity]) -> List[Cavity]:
        """Create a list of Cavity instances from a list of Mole2Cavity instances.

        Args:
            stru (Structure): The structure instance containing the cavity.
            mole2cavities (List[Mole2Cavity]): A list of Mole2Cavity Instance to initialize Cavity instances.
        """
        cavities = [cls.from_mole2_cavity(stru, mole2cavity) for mole2cavity in mole2cavities]
        return cavities

    
    def __eq__(self, other: Cavity) -> bool:
        return self.stru == other.stru and set(self.boundary_residues) == set(other.boundary_residues) and set(self.inner_residues) == set(other.inner_residues)

def identify_stru_cavities(stru: Structure,
        non_active_residues: List[Residue] = list(), 
        probe: float = None, 
        inner: float = None, 
        mesh_density: float = None,
        ignore_hetatm: bool = None,
        work_dir: str = None,
        use_mono: bool = True
    ) -> List[Cavity]:
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
    sp.save_structure(outfile=pdb_filepath, stru=stru)
    non_active_parts = [resi.key() for resi in non_active_residues]

    mole2cavities = mole2.identify_cavities(pdb_path=pdb_filepath, non_active_parts=non_active_parts,
        probe=probe, inner=inner, mesh_density=mesh_density, ignore_hetatm=ignore_hetatm,
        work_dir=work_dir, use_mono=use_mono)
    
    stru_cavities = Cavity.from_mole2_cavities(stru=stru, mole2cavities=mole2cavities)
    fs.safe_rm(pdb_filepath)
    return stru_cavities

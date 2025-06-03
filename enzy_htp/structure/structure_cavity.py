"""This class defines the cavity of a Structure instance.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>

Date: 2025-06-02
"""

from __future__ import annotations
from typing import List, Tuple

import numpy as np
import numpy.typing as npt
import pyvista as pv

from enzy_htp import interface, _LOGGER, PDBParser, Structure
from enzy_htp._interface import Mole2Cavity

class StructureCavity(Mole2Cavity):
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
        self.boundary_residues = [stru.find_residue_with_key(key) for key in self.boundary_residue_keys]
        self.inner_residues = [stru.find_residue_with_key(key) for key in self.inner_residue_keys]
        return
    
    @classmethod
    def from_mole2_cavity(cls, stru: Structure, mole2cavity: Mole2Cavity) -> StructureCavity:
        """Create a Cavity instance from a Mole2Cavity instance.
        
        Args:
            stru (Structure): The structure instance containing the cavity.
            mole2cavity (Mole2Cavity): The Mole2Cavity Instance to initialize the Cavity instance.
        """
        cavity = cls.__init__(stru, points=mole2cavity.points(), 
                mesh=mole2cavity.mesh_, probe=mole2cavity.probe(), inner=mole2cavity.inner(), 
                mesh_density=mole2cavity.mesh_density(), com=mole2cavity.com_, mole2_volume=mole2cavity.mole2_volume, 
                boundary_residue_keys=mole2cavity.boundary_residue_keys, inner_residue_keys=mole2cavity.inner_residue_keys)
        return cavity
    
    @classmethod
    def from_mole2_cavities(cls, stru: Structure, mole2cavities: List[Mole2Cavity]) -> List[StructureCavity]:
        """Create a list of Cavity instances from a list of Mole2Cavity instances.

        Args:
            stru (Structure): The structure instance containing the cavity.
            mole2cavities (List[Mole2Cavity]): A list of Mole2Cavity Instance to initialize Cavity instances.
        """
        cavities = [cls.from_mole2_cavity(stru, mole2cavity) for mole2cavity in mole2cavities]
        return cavities

    
    def __eq__(self, other: StructureCavity) -> bool:
        return self.stru == other.stru and set(self.boundary_residues) == set(other.boundary_residues) and set(self.inner_residues) == set(other.inner_residues)

#! python3
# -*- coding: utf-8 -*-
"""
@File   : cavity.py
@Created: 2025/06/22 17:41
@Author : Zhong, Yinjie
@Email  : yinjie.zhong@vanderbilt.edu
"""

from __future__ import annotations
from typing import List, Tuple, Literal
from .. import Structure, Residue
import pyvista as pv
import numpy.typing as npt

class Cavity():
    """Class for the cavity of a Structure that stores information about individual cavity identified. 
    Supported operations include cavity volume and center-of-mass calculations, 
    and indicating if the 3D cartesian point is contained within the cavity.
    """
    def __init__(self, stru: Structure, mesh: pv.PolyData, probe: float, 
            inner: float, mesh_density: float, software_report_volume: float,            
            boundary_residues: List[Residue], inner_residues: List[Residue],
            type: Literal["Cavity", "Void", None] = None):
        """Initialize a Cavity instance.
        
        Args:
            stru (Structure): The structure instance containing the cavity.
            mesh (pv.PolyData): The mesh representing the cavity established with pyvista.
            probe (float): Probe radius used during collection in Angstroms.
            inner (float): Inner radius used during collection in Angstroms.
            software_report_volume (float): The volume of the cavity in A^3 calculated by software or engine.
            boundary_residues (List[Residue]): A list() of residues forming the cavity that are on the boundary of the structure.
            inner_residues (List[Residue]): A list() of residues forming the cavity that are inside the structure.
            type (str, optional): The type of the cavity.
                * "Cavity": This cavity connects to the surface of the structure;
                * "Void": This cavity is hidden inside the structure.
                * None: The type value is undefined or not provided.
        """
        self.stru = stru
        self.mesh = mesh
        self.probe = probe
        self.inner = inner
        self.mesh_density = mesh_density
        self.software_report_volume = software_report_volume
        self.boundary_residues = boundary_residues
        self.inner_residues = inner_residues
        self.type = type
        return
    
    @property
    def points(self):
        """The points of the mesh as a numpy array with format (n, 3)."""
        return self.mesh.points
    
    @property
    def volume(self):
        """The volume of the cavity calculated from the mesh."""
        return self.mesh.volume
    
    @property
    def center_of_mass(self):
        """The center-of-mass of the mesh as a numpy array with format (x, y, z)."""
        return self.mesh.center_of_mass()

    def contains(self, point: npt.NDArray) -> bool:
        """Check if the given point is inside the cavity mesh (turn off surface checking).

        Args:
            point (npt.NDArray): The point in question as a numpy array with format (x, y, z).

        Returns:
            If the point is contained by the mesh.
        """
        points_poly = pv.PolyData([point, ])
        result = points_poly.select_enclosed_points(
            self.mesh,
            check_surface=False
        )
        return bool(result['SelectedPoints'][0])
    
    def __eq__(self, other: Cavity) -> bool:
        return self.stru == other.stru and set(self.boundary_residues) == set(other.boundary_residues) and set(self.inner_residues) == set(other.inner_residues)

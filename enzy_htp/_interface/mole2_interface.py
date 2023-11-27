"""Defines a Mole2Interface class that serves as a bridge for enzy_htp to utilize Mole2. Uses the Mole2Config
class found in enzy_htp/_config/mole2_config.py. Supported operations include:

    + getting cavities from a given structure

Sub-module includes a companion Mole2Cavity class that represents a cavity in a structure. The Mole2Cavity objects
stores information about the cavity. Supported operations include:

    + calculation of cavity volume
    + indication if 3D cartesian points are included in a cavity.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-26
"""

from pathlib import Path
from typing import List, Tuple

import numpy as np
import numpy.typing as npt
import pyvista as pv

from .base_interface import BaseInterface

from enzy_htp import config as eh_config
from enzy_htp import _LOGGER 
from enzy_htp.core import file_system as fs

from enzy_htp._config.mole2_config import Mole2Config, default_mole2_config

class Mole2Cavity:
    """Companion class to the Mole2Interface that stores information about individual cavities identified
    by Mole2. Supported operations include cavity volume and center of mass calculations, and indicating if
    3D cartesian points are contained within the cavity.

    Attributes:
        points_ : A list() of 3D vertices describing the cavity.
        mesh_ : A pyvista.PolyData object representing the cavity's mesh.
        probe_ : Probe radius used during collection in Angstroms.
        inner_ : Inner radius used during collection in Angstroms.
        mesh_density_ : Mesh density used during collection in Angstroms.
        com_ : The center-of-mass of the mesh as a numpy array with format (x, y, z).
    """ 

    def __init__(self, points :List[npt.NDArray], mesh : pv.PolyData, probe:float, inner:float, mesh_density:float, com:npt.NDArray):
        """Simplistic constructor. Each attribute is directly set by a parameter."""
        self.points_ = points
        self.mesh_ = mesh
        self.probe_ = probe
        self.inner_ = inner
        self.mesh_density_ = mesh_density
        self.com_ = com 

    
    def points(self) -> List[npt.NDArray]:
        """TODO(CJ)"""
        return self.points_

    def volume(self) -> float:
        """The volume of the cavity in A^3. Set at time of construction."""
        return self.mesh_.volume

    def center_of_mass(self) -> npt.NDArray:
        """The center-of-mass of the mesh as a numpy array with format (x, y, z). Set at time of construction."""
        return self.com_

    def contains(self, point: npt.NDArray) -> bool:
        """Does the mesh contain a given point? Uses pyvista.PolyData and turns off surface checking.

        Args:
            point: The point in question as a numpy array with format (x, y, z).

        Returns:
            If the point is contained by the mesh.
        """

        point = pv.PolyData([point])
        result = point.select_enclosed_points(
            self.mesh_,
            check_surface=False
        )
        
        return bool(result['SelectedPoints'][0])

    def contains_points(self, points: List[npt.NDArray]) -> bool:
        """TODO(CJ): make this a dispatch and this will be the list version"""

        points = pv.PolyData(points)
        result = points.select_enclosed_points(
            self.mesh_,
            check_surface=False
        )
        
        return result['SelectedPoints']


    def probe(self) -> float:
        """Getter for the probe radius in A"""
        return self.probe_

    def inner(self) -> float:
        """Getter for the inner radius in A"""
        return self.inner_
    
    def mesh_density(self) -> float:
        """Getter for the mesh density in A"""
        return self.mesh_density_

class Mole2Interface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize Mole2. Supported operations
    include filling structures with ligand transplants. Users should use this class as the only way to 
    interact with this application.

    Attributes:
        config_ : The Mole2Config() class which provides settings for both running Mole2 and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """
    def __init__(self, parent, config: Mole2Config = None) -> None:
        """Simplistic constructor that optionally takes a Mole2Config object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_mole2_config)


    def _write_xml_input(self,  molfile:str,
                                work_dir:str,
                                non_active_parts:List[Tuple[str,int]], 
                                probe:float, 
                                inner:float, 
                                mesh_density:float,
                                ignore_hetatm:bool
                                ) -> str:
        """Given settings for the mole2 run, writes an input .xml file for mole2. Always writes the file 
        to work_dir/__mole2_input.xml.

        Args:
            molfile: The .pdb file to look for cavities in.
            work_dir: Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
            non_active_parts: Residues that should be skipped. Should be in format List[Tuple] where Tuple has format (chain id, residue number).
            probe: Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
            inner: Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
            mesh_density: Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
            ignore_hetam: TOOD(CJ)

        Returns:
            A str() with the filename of the .xml containing the mole2 input.
        """
        content:List[str] = [
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<Tunnels>",
           f"\t<WorkingDirectory>{work_dir}/</WorkingDirectory>",
           f"\t<Input>{molfile}</Input>\n"]

        if non_active_parts is not None:
            content.append("\t<NonActiveParts>")
            for (chain, rnum) in non_active_parts:
                content.append(f"\t\t<Residue Chain=\"{chain}\" SequenceNumber=\"{rnum}\" />")
        
            content.append("\t</NonActiveParts>")

        content.extend([
            "\t<Params>",
           f"\t\t<Cavity ProbeRadius=\"{probe}\" InteriorThreshold=\"{inner}\" IgnoreHETAtoms=\"{ignore_hetatm}\"/>",
            "\t</Params>",
            "\t<Export>",
            "\t\t<Formats Mesh=\"1\" />",
            "\t\t<Mesh Density=\"0.5\" />",
            "\t</Export>",
            "</Tunnels>",
        ])

        outfile:str = f"{work_dir}/__mole2_input.xml"

        fs.write_lines(outfile, content)

        return outfile

    def _parse_cavity(self, fname:str, probe:float, inner:float, mesh_density:float) -> Mole2Cavity:
        """Factory function to produce Mole2Cavity objects. Each object stores information about both the
        cavity itself and the settings used to collect it.

        Args:
            fname: The .cavity file from a mole2 run.
            probe: Probe radius used during collection in A.
            inner: Inner radius used during collection in A.
            mesh_density: Mesh density used during collection in A.
        
        Returns:
            A newly constructed Mole2Cavity.
        """
        
        lines:List[str] = fs.lines_from_file(fname)
        points = []
        num_lines = int(lines.pop(0))
        
        for ll in lines[:num_lines]:
            points.append(np.array(ll.split()).astype(float))
        
        lines = lines[num_lines:]
        
        points = np.array( points )
        num_pgons = int(lines.pop(0))
        
        
        cnct = list(map(int,lines))
        cnct.reverse()
        surface = []
        
        while cnct:
            n = cnct.pop()
        
            surface.append( [n] + [cnct.pop() for _ in range(n)] )
        com = 0.0
        verts = list()

        for ss in surface:
            for idx in ss[1:-1]:
                verts.append( points[idx] )
        
        verts = np.array(verts)
        com = np.mean(verts,axis=0)

        return Mole2Cavity(
            points=points,
            mesh=pv.PolyData(points, np.hstack(surface)),
            probe=probe,
            inner=inner,
            mesh_density=mesh_density,
            com=com
            )

    def identify_cavities(self, molfile:str, 
                                non_active_parts:List[Tuple[str,int]]=None, 
                                probe:float=None, 
                                inner:float=None, 
                                mesh_density:float=None,
                                ignore_hetatm:bool=None,
                                work_dir:str=None,
                                use_mono:bool=True
                                ) -> List[Mole2Cavity]:
        """Identifies cavities in a protein structure using the Mole2 software package. Client method that should be 
        called by users. Results are represented via Mole2Cavity objects that support basic geometry operations.

        Args:
            molfile: The .pdb file to look for cavities in.
            non_active_parts: Residues that should be skipped. Should be in format List[Tuple] where Tuple has format (chain id, residue number).
            probe: Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
            inner: Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
            mesh_density: Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
            ignore_hetatm: TODO(CJ)
            work_dir: Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
            use_mono: Does mono need to be used during run time? Defaults to true.

        Returns:
            A list() of Mole2Cavity ebjects.
        """

        if probe is None:
            probe = self.config_.PROBE

        if inner is None:
            inner = self.config_.INNER

        if mesh_density is None:
            mesh_density = self.config_.MESH_DENSITY

        if work_dir is None:
            work_dir = eh_config['system.SCRATCH_DIR']

        if ignore_hetatm is None:
            ignore_hetatm = self.config_.IGNORE_HETATM

        fs.safe_mkdir( work_dir )

        fs.safe_rmdir( f"{work_dir}/mesh/")

        fs.check_not_empty( molfile )

        xml_file:str = self._write_xml_input(molfile, work_dir, non_active_parts, probe, inner, mesh_density, ignore_hetatm)

        if use_mono:
            self.env_manager_.run_command(self.config_.MONO, [self.config_.MOLE2, xml_file])
        else:
            self.env_manager_.run_command(self.config_.MOLE2, [xml_file])
        
        mesh_files:List[str] = list(Path(f"{work_dir}/mesh/").glob("cavity_*.mesh"))
        _LOGGER.info(f"Found {len(mesh_files)} cavities using probe radius of {probe:.3f} A and inner radius of {inner:.3f} A")
            
        result:List[Mole2Cavity] = list()
        for mf in mesh_files:
            result.append(self._parse_cavity(mf, probe, inner, mesh_density))
                    
        return result

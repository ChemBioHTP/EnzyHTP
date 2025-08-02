"""Defines a Mole2Interface class that serves as a bridge for enzy_htp to utilize Mole2. Uses the Mole2Config
class found in enzy_htp/_config/mole2_config.py. Supported operations include:

    + getting cavities from a given structure

A cavity in a structure is represented by [Cavity class](../structure/structure_cavity/cavity.py).
The Mole2Cavity objects stores information about the cavity. Supported operations include:

    + calculation of cavity volume
    + indication if 3D cartesian points are included in a cavity.

Author:
    - Chris Jurich <chris.jurich@vanderbilt.edu>;
    - Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>;
Date: 2025-06-24
"""
from os import path
from pathlib import Path
from typing import List, Tuple, Literal
import xml.etree.ElementTree as ET

import numpy as np
import numpy.typing as npt
import pyvista as pv

from enzy_htp.structure import Structure, Residue, Cavity

from .base_interface import BaseInterface

from enzy_htp import _LOGGER, config as eh_config, PDBParser
from enzy_htp.core import file_system as fs
from enzy_htp.structure import Structure, Residue, Cavity

from enzy_htp._config.mole2_config import Mole2Config, default_mole2_config

sp = PDBParser()

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


    def _write_xml_input(self,  pdb_path: str,
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
            pdb_path: The .pdb file to look for cavities in.
            work_dir: Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
            non_active_parts: Residues that should be skipped. Should be in format List[Tuple] where Tuple has format (chain id, residue number).
            probe: Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
            inner: Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
            mesh_density: Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
            ignore_hetam: TOOD(CJ)

        Returns:
            input_filepath (str): The path to the .xml file containing the mole2 input.
        """
        content:List[str] = [
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<Tunnels>",
           f"\t<WorkingDirectory>{work_dir}/</WorkingDirectory>",
           f"\t<Input>{pdb_path}</Input>\n"]

        if non_active_parts is not None:
            content.append("\t<NonActiveParts>")    # The `NonActiveResidues` described in the official doc is not correct.
            for (chain, rnum) in non_active_parts:
                content.append(f"\t\t<Residue Chain=\"{chain}\" SequenceNumber=\"{rnum}\" />")
        
            content.append("\t</NonActiveParts>")

        content.extend([
           f"\t<Params>",
           f"\t\t<Cavity IgnoreHETAtoms=\"{ignore_hetatm}\" ProbeRadius=\"{probe}\" InteriorThreshold=\"{inner}\"/>",
            "\t</Params>",
            "\t<Export>",
            "\t\t<Formats Mesh=\"1\" />",
            "\t\t<Mesh Density=\"1.33\" />",
            "\t</Export>",
            "</Tunnels>",
        ])

        input_filepath: str = path.join(work_dir, "mole2_input.xml")

        fs.write_lines(input_filepath, content)

        return input_filepath

    def _read_cavity_from_xml(self, cavities_xml_filepath: str, cavity_id: int, cavity_type: Literal["Cavity", "Void"] = "Cavity") -> Tuple[float, list, list]:
        """Given a .xml file from a mole2 run, parses the cavity information and returns the boundary
        and inner residues information.
        
        Args:
            cavities_xml_filepath (str): The .xml filepath from a mole2 run.
            cavity_id (int): The cavity id.
            
        Returns:
            A tuple containing (volume, boundary_residue_keys, inner_residue_keys).
        """
        if (not cavities_xml_filepath or not cavity_id):
            return None, None, None
        tree = ET.parse(cavities_xml_filepath)
        root = tree.getroot()
        
        # Find the cavity with matching Id
        cavity = root.find(f".//Cavity[@Type='{cavity_type}'][@Id='{cavity_id}']")
        if cavity is None:
            raise ValueError(f"Cavity with Id={cavity_id} not found in {cavities_xml_filepath}")
        
        mole2_volume = None
        try:
            mole2_volume = float(cavity.attrib.get("Volume"))
        except ValueError as exc:
            raise ValueError(f"Could not parse cavity volume from {cavities_xml_filepath}") from exc
        
        # Get Boundary and Inner residues text
        boundary_residues_text = cavity.find(".//Boundary/Residues").text or ""
        boundary_residues_text_split = boundary_residues_text.split(",") if boundary_residues_text else list()
        inner_residues_text = cavity.find(".//Inner/Residues").text or ""
        inner_residues_text_split = inner_residues_text.split(",") if inner_residues_text else list()

        # Parse the text into a list of residues
        boundary_residue_keys = [(txt.split()[-1], int(txt.split()[-2])) for txt in boundary_residues_text_split]
        inner_residue_keys = [(txt.split()[-1], int(txt.split()[-2])) for txt in inner_residues_text_split]
        
        return mole2_volume, boundary_residue_keys, inner_residue_keys

    def _parse_cavity(self, stru: Structure, mesh_filepath: str, probe: float, inner: float, 
            mesh_density: float, cavity_id: int = None, 
            cavity_xml_filepath: str = None, cavity_type: Literal["Cavity", "Void"] = "Cavity") -> Cavity:
        """Factory function to produce Mole2Cavity objects. Each object stores information about both the
        cavity itself and the settings used to collect it.

        Args:
            stru (Structure): The structure instance to detect cavities from.
            mesh_filepath: The .mesh file from a mole2 run.
            probe: Probe radius used during collection in A.
            inner: Inner radius used during collection in A.
            mesh_density: Mesh density used during collection in A.
        
        Returns:
            A newly constructed Mole2Cavity.
        """
        
        lines:List[str] = fs.lines_from_file(mesh_filepath) # Read cavity_X.mesh file lines.
        points = []
        num_lines = int(lines.pop(0))       # The first line is the number of points.
        
        for ll in lines[:num_lines]:        # Read point coordinates by number.
            points.append(np.array(ll.split()).astype(float))
        points = np.array(points)           # Convert to numpy array
        
        lines = lines[num_lines:]       # The remaining rows are surface data.
        num_pgons = int(lines.pop(0))   # Surface number.
        
        cnct = list(map(int,lines))     # Convert face data to integer list
        cnct.reverse()                  # Reverse for later processing

        surfaces = []
        while cnct:
            n = cnct.pop()              # Get number of vertices in surface.
            surfaces.append([n-1] + [cnct.pop() for _ in range(n)][:-1])    # Store [number of vertices + vertex index]
            continue
    
        com = 0.0
        verts = list()
        for ss in surfaces:
            for idx in ss[1:-1]:            # Get the vertex index of the face (skip the first element: the vertex number).
                verts.append(points[idx])   # Collect vertex coordinates.
        
        verts = np.array(verts)
        com = np.mean(verts,axis=0)     # Calculate the geometric center of the cavity.

        # print(f"Surface elements: {len(surfaces)}, Surface count: {num_pgons}")
        mesh = pv.PolyData(var_inp=points, faces=np.hstack(surfaces))
        mesh = mesh.clean()
        mesh = mesh.triangulate()

        mole2_volume, boundary_residue_keys, inner_residue_keys = self._read_cavity_from_xml(cavity_xml_filepath, cavity_id, cavity_type)
        boundary_residues = [stru.find_residue_with_key(key) for key in boundary_residue_keys] if boundary_residue_keys else []
        inner_residues = [stru.find_residue_with_key(key) for key in inner_residue_keys] if inner_residue_keys else []
        return Cavity(
            stru=stru,
            mesh=mesh,          # Building a grid with pyvista
            probe=probe,
            inner=inner,
            mesh_density=mesh_density,
            software_report_volume=mole2_volume,
            boundary_residues=boundary_residues,
            inner_residues=inner_residues,
            cavity_type=cavity_type
        )

    def identify_cavities(self, stru: Structure, 
            non_active_residues: List[Residue] = list(), 
            probe: float = None, 
            inner: float = None, 
            mesh_density: float = None,
            ignore_hetatm: bool = None,
            work_dir: str = None,
            use_mono: bool = True,
            **kwargs
        ) -> List[Cavity]:
        """Identifies cavities in a protein structure using the Mole2 software package. Client method that should be 
        called by users. Results are represented via Mole2Cavity objects that support basic geometry operations.

        Args:
            stru (Structure): The structure instance to detect cavities from.
            non_active_residues (List[Residue], optional): Residues that should be skipped.
            probe (float, optional): Probe radius to use in A. Defaults to Mole2Config.PROBE if not supplied.
            inner (float, optional): Inner radius to use in A. Defaults to Mole2Config.INNER if not supplied.
            mesh_density (float, optional): Mesh density to use in A. Defaults to Mole2Config.MESH_DENSITY if not supplied.
            ignore_hetatm (bool, optional): TODO (CJ)
            work_dir (str, optional): Directory to do work in. Defaults to system.SCATCH_DIR if not supplied.
            use_mono (bool, optional): Indicate if mono need to be used during run time. Defaults to true.

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

        non_active_parts = [resi.key() for resi in non_active_residues]

        fs.safe_mkdir(work_dir)
        fs.safe_rmdir(f"{work_dir}/mesh/")
        
        pdb_filepath = fs.get_valid_temp_name(path.join(work_dir, "stru_cavity_temp.pdb"), ext_set=["pdb"])
        sp.save_structure(outfile=pdb_filepath, stru=stru)

        input_xml_filepath: str = self._write_xml_input(pdb_filepath, work_dir, non_active_parts, probe, inner, mesh_density, ignore_hetatm)

        if use_mono:
            self.env_manager_.run_command(self.config_.MONO, [self.config_.MOLE2, input_xml_filepath])
        else:
            self.env_manager_.run_command(self.config_.MOLE2, [input_xml_filepath])
        
        cavity_mesh_files: List[str] = [str(filepath.resolve()) for filepath in Path(f"{work_dir}/mesh/").glob("cavity_*.mesh")]
        void_mesh_files: List[str] = [str(filepath.resolve()) for filepath in Path(f"{work_dir}/mesh/").glob("void_*.mesh")]
        cavities_xml_file = str(Path(work_dir).joinpath("xml", "cavities.xml"))
        _LOGGER.info(f"Found {len(cavity_mesh_files)} cavities and {len(void_mesh_files)} void cavities using probe radius of {probe:.3f} A and inner radius of {inner:.3f} A")
        
        result: List[Cavity] = list()
        for i, mf in enumerate(cavity_mesh_files):
            result.append(
                self._parse_cavity(stru=stru, mesh_filepath=mf, probe=probe, inner=inner, mesh_density=mesh_density, 
                    cavity_id=(i+1), cavity_xml_filepath=cavities_xml_file, cavity_type="Cavity")
            )
        for i, mf in enumerate(void_mesh_files):
            result.append(
                self._parse_cavity(stru=stru, mesh_filepath=mf, probe=probe, inner=inner, mesh_density=mesh_density, 
                    cavity_id=(i+1), cavity_xml_filepath=cavities_xml_file, cavity_type="Void")
            )

        fs.clean_temp_file_n_dir([cavity_mesh_files, void_mesh_files, cavities_xml_file, input_xml_filepath])
        
        return result

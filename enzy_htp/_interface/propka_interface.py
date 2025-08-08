"""Submodule contains code for interfacing with the PROPKA software package.
+ get_residue_pka_from_pdb
    Run PROPKA on a pdb file and parse the output to get the pKa of each residue.
+ get_residue_pka_from_stru
    Run PROPKA on a Structure object and return the pKa values.

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-07-18
"""
from typing import Dict, List, Any, Tuple, Union
import os
import tempfile
import importlib

from .base_interface import BaseInterface
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.file_system import safe_mkdir
from enzy_htp._config.propka_config import PropkaConfig, default_propka_config
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp import config as eh_config

# Importing PROPKA components
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer


class PropkaInterface(BaseInterface):
    """Interface class for PROPKA software package."""

    def __init__(self, parent, config: PropkaConfig = None):
        """Constructor for PropkaInterface."""
        super().__init__(parent, config, default_propka_config)
    
    def set_parent(self, parent):
        """Set the parent interface."""
        self.parent_ = parent

    def missing_executables(self) -> List[str]:
        return []

    def missing_env_vars(self) -> List[str]:
        return []

    def missing_py_modules(self) -> List[str]:
        missing = []
        for module in ['propka']:
            if importlib.util.find_spec(module) is None:
                missing.append(module)
        return missing

    def get_residue_pka_from_stru(self, stru: Structure, work_dir: str = None) -> Dict[Tuple[str, int], float]:
        """Calculate pKa values for residues in a Structure object using PROPKA.

        Args:
            stru: The Structure object to analyze.
            work_dir: The working directory for the calculation. If None, uses system SCRATCH directory.

        Returns:
            Dict[Tuple[str, int], float]: A dictionary mapping (chain_id, residue_number) 
                                         tuples to their pKa values.
        """
        # Set default work_dir to SCRATCH if not provided
        if work_dir is None:
            work_dir = eh_config.system.SCRATCH_DIR
            
        safe_mkdir(work_dir)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False, dir=work_dir) as tmp_file:
            temp_pdb_path = tmp_file.name
            PDBParser().save_structure(temp_pdb_path, stru)
        
        try:
            detailed_results = self.run_propka(temp_pdb_path, work_dir)
            # Convert detailed results to simple dict for backward compatibility
            simple_results = {}
            for res_data in detailed_results:
                key = (res_data["chain_id"], res_data["res_num"])
                simple_results[key] = res_data["pKa"]
            return simple_results
        finally:
            if os.path.exists(temp_pdb_path):
                os.unlink(temp_pdb_path)

    def run_propka(self, pdb_path: str, work_dir: str = None) -> List[Dict[str, Any]]:
        """
        Main wrapper for PROPKA that calculates pKa values for residues in a given PDB file.

        Args:
            pdb_path (str): The path to the PDB file.
            work_dir (str, optional): The working directory for the calculation. If None, uses system SCRATCH directory.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing detailed pKa information 
                                 for each ionizable residue.
        """
        try:
            # Set default work_dir to SCRATCH if not provided
            if work_dir is None:
                work_dir = eh_config.system.SCRATCH_DIR
                
            safe_mkdir(work_dir)
            original_cwd = os.getcwd()
            os.chdir(work_dir)
            
            # PROPKA calculation
            options = loadOptions([os.path.abspath(pdb_path)])
            parameters = read_parameter_file(options.parameters, Parameters())
            my_molecule = MolecularContainer(parameters, options)
            my_molecule = read_molecule_file(os.path.abspath(pdb_path), my_molecule)
            my_molecule.calculate_pka()

            # Extract detailed pKa information
            conformation = my_molecule.conformations.get('AVR', {})
            residues_data = []
            for group in conformation.groups:
                atom = group.atom
                row_dict = {}
                row_dict["res_num"] = atom.res_num
                row_dict["ins_code"] = atom.icode
                row_dict["res_name"] = atom.res_name
                row_dict["chain_id"] = atom.chain_id
                row_dict["group_label"] = group.label
                row_dict["group_type"] = getattr(group, "type", None)
                row_dict["pKa"] = group.pka_value
                row_dict["model_pKa"] = group.model_pka
                row_dict["buried"] = group.buried
                if group.coupled_titrating_group:
                    row_dict["coupled_group"] = group.coupled_titrating_group.label
                else:
                    row_dict["coupled_group"] = None
                residues_data.append(row_dict)

            return residues_data

        except FileNotFoundError:
            _LOGGER.error(f"PDB file not found: {pdb_path}")
            raise
        except Exception as e:
            _LOGGER.error(f"Error in get_residue_pka: {e}")
            raise
        finally:
            os.chdir(original_cwd)


propka_interface = PropkaInterface(None, eh_config._propka)
"""The singleton of PropkaInterface() that handles all PROPKA related operations in EnzyHTP.
Instantiated here so that other _interface subpackages can use it."""
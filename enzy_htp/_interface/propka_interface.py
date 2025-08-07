"""Submodule contains code for interfacing with the PROPKA software package.
+ get_residue_pka_from_pdb
    Run PROPKA on a pdb file and parse the output to get the pKa of each residue.
+ get_residue_pka_from_stru
    Run PROPKA on a Structure object and return the pKa values.

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-07-18
"""
from typing import Dict, List
import os
import tempfile
import importlib

from .base_interface import BaseInterface
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.file_system import safe_mkdir
from enzy_htp._config.propka_config import PropkaConfig, default_propka_config
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
# remove_solvent imported dynamically to avoid circular import
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

    def get_residue_pka_from_stru(self, stru: Structure, work_dir: str = "./propka", remove_solvents: bool = True) -> Dict[int, float]:
        """Calculate pKa values for residues in a Structure object using PROPKA.

        Args:
            stru: The Structure object to analyze.
            work_dir: The working directory for the calculation.
            remove_solvents: Whether to remove solvent molecules before calculation.
                           Default True since PROPKA doesn't benefit from solvents.

        Returns:
            Dict[int, float]: A dictionary mapping residue numbers to their pKa values.
        """
        safe_mkdir(work_dir)
        
        # Remove solvents by default as PROPKA doesn't benefit from them
        if remove_solvents:
            from enzy_htp.preparation.clean import remove_solvent
            stru_clean = remove_solvent(stru, in_place=False)
        else:
            stru_clean = stru
        
        # Use omit_chain_id=True if we still have many residues after solvent removal
        # to prevent chain ID overflow issues
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False, dir=work_dir) as tmp_file:
            temp_pdb_path = tmp_file.name
            PDBParser.save_structure(temp_pdb_path, stru_clean, omit_chain_id=True)
        
        try:
            result = self.get_residue_pka_from_pdb(temp_pdb_path, work_dir)
            return result
        finally:
            if os.path.exists(temp_pdb_path):
                os.unlink(temp_pdb_path)

    def get_residue_pka_from_pdb(self, pdb_path: str, work_dir: str = "./propka") -> Dict[int, float]:
        """
        Calculates the pKa values for residues in a given PDB file using PROPKA.

        Args:
            pdb_path (str): The path to the PDB file.
            work_dir (str, optional): The working directory for the calculation. Defaults to "./propka".

        Returns:
            Dict[int, float]: A dictionary mapping residue numbers to their pKa values.
        """
        try:
            safe_mkdir(work_dir)
            original_cwd = os.getcwd()
            os.chdir(work_dir)
            
            # PROPKA calculation
            options = loadOptions([os.path.abspath(pdb_path)])
            parameters = read_parameter_file(options.parameters, Parameters())
            my_molecule = MolecularContainer(parameters, options)
            my_molecule = read_molecule_file(os.path.abspath(pdb_path), my_molecule)
            my_molecule.calculate_pka()

            # Extract pKa values
            conformation = my_molecule.conformations.get('AVR', {})
            residues_pka = {}
            for group in conformation.groups:
                atom = group.atom
                residues_pka[atom.res_num] = group.pka_value

            return residues_pka

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
"""Submodule contains code for interfacing with the PROPKA software package.
+ get_residue_pka_from_pdb
    Run PROPKA on a pdb file and parse the output to get the pKa of each residue.

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-07-18
"""
from typing import Dict, List
import os
import importlib
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.file_system import safe_mkdir

# Importing PROPKA components
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer

def missing_executables() -> List[str]:
    return []

def missing_env_vars() -> List[str]:
    return []

def missing_py_modules() -> List[str]:
    missing = []
    for module in ['propka']:
        if importlib.util.find_spec(module) is None:
            missing.append(module)
    return missing

def get_residue_pka_from_stru(stru, work_dir: str = "./propka") -> Dict[int, float]:
    pass

def get_residue_pka_from_pdb(pdb_path: str, work_dir: str = "./propka") -> Dict[int, float]:
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

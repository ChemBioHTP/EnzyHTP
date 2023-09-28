"""Defines PyMolConfig() which holds configuration settings for enzy_htp to interface with the 
PyMol software package. 
These configurations are project level configurations, mainly composed by interfacing related configs,
default settings of specific PyMol functions, default resource settings of PyMol related jobs, 
and non-task specific settings like output level etc. 
(NOTE: decouple settings actually used in specific task with configs here) 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2023-02-16
"""
from pprint import pprint
from copy import deepcopy
from typing import Any, List, Dict, Tuple

from .base_config import BaseConfig


class PyMolConfig(BaseConfig):
    """Class that holds default values for running PyMol within enzy_htp."""

    WGET: str = "wget"
    """ """

    PYMOL_TO_ATOM_MAPPER: Dict[str, str] = {
        "1HB": "HB1",
        "2HB": "HB2",
        "3HB": "HB3",
        "1HG": "HG1",
        "2HG": "HG2",
        "3HG": "HG3",
        "1HD": "HD1",
        "2HD": "HD2",
        "3HD": "HD3",
        "1HE": "HE1",
        "2HE": "HE2",
        "3HE": "HE3",
        "2HH": "HH2",
        "1HZ": "HZ1",
        "2HZ": "HZ2",
        "3HZ": "HZ3",
        "1HH1": "HH11",
        "2HH1": "HH12",
        "1HH2": "HH21",
        "2HH2": "HH22",
        "1HE2": "HE21",
        "2HE2": "HE22",
        "1HG2": "HG21",
        "2HG2": "HG22",
        "3HG2": "HG23",
        "1HG1": "HG11",
        "2HG1": "HG12",
        "3HG1": "HG13",
        "1HD1": "HD11",
        "2HD1": "HD12",
        "3HD1": "HD13",
        "1HD2": "HD21",
        "2HD2": "HD22",
        "3HD2": "HD23",
    }
    """A map to convert messed up PyMOL atom names to a more conventional name for the Structure class."""


    DEFAULT_OUTPUT_LV: Tuple[str] = ("disable", "all", "everything")
    """the default output level of a new pymol session in enzy_htp: mute everything.
    used in session.cmd.feedback() ref: https://pymolwiki.org/index.php/Feedback"""

    STRUCTURE_STEM: str = "https://files.rcsb.org/download/"
    """ """

    LIGAND_STEM: str = "https://files.rcsb.org/ligands/download/"
    """ """

    def display(self) -> None:
        """TODO"""
        pass

    def required_env_vars(self) -> List[str]:
        """ """
        return list()

    def required_executables(self) -> List[str]:
        """ """
        return [self.WGET]

    def required_py_modules(self) -> List[str]:
        """ """
        return ["pymol2"]
    
    def get_canonical_atom_name(self, res: str, atom: str) -> str:
        """Updates PyMOL atom name to name consistent with Structure class.
        Args:
            res: the current residue.
            atom: the current atom.
        Returns:
            The correct atom name. 
        """
        if res == "GLY" and atom == "HA":
            return "HA2"

        res = self.PYMOL_TO_ATOM_MAPPER.get(atom, None)
        return res

def default_pymol_config() -> PyMolConfig:
    """Creates a deep-copied default version of the PyMolConfig() class."""
    return deepcopy(PyMolConfig())

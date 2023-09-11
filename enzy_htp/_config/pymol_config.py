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

    WGET:str="wget"
    """ """

    DEFAULT_OUTPUT_LV: Tuple[str] = ("disable", "all", "everything")
    """the default output level of a new pymol session in enzy_htp: mute everything.
    used in session.cmd.feedback() ref: https://pymolwiki.org/index.php/Feedback"""

    STRUCTURE_STEM:str="https://files.rcsb.org/download/"
    """ """

    LIGAND_STEM:str="https://files.rcsb.org/ligands/download/"
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


def default_pymol_config() -> PyMolConfig:
    """Creates a deep-copied default version of the PyMolConfig() class."""
    return deepcopy(PyMolConfig())

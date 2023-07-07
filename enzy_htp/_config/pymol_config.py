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


class PyMolConfig:
    """Class that holds default values for running PyMol within enzy_htp."""

    DEFAULT_OUTPUT_LV: Tuple[str] = ("disable", "all", "everything")
    """the default output level of a new pymol session in enzy_htp: mute everything.
    used in session.cmd.feedback() ref: https://pymolwiki.org/index.php/Feedback"""

    def __init__(self, parent=None):
        """Trivial constructor that optionally sets parent_ dependency. parent_ is None by default."""
        self.parent_ = parent

    def display(self) -> None:
        """TODO"""
        pass

    def __getitem__(self, key: str) -> Any:
        """Getter that enables [] accession of PyMolConfig() attributes."""
        if key.count("."):
            key1, key2 = key.split(".", 1)
            return getattr(self, key1)[key2]
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter that enables [] accession of PyMolConfig() attributes with value validation."""
        if key.count("."):
            key1, key2 = key.split(".")
            PyMolConfig.__dict__[key1][key2] = value
        else:
            setattr(self, key, value)


def default_pymol_config() -> PyMolConfig:
    """Creates a deep-copied default version of the PyMolConfig() class."""
    return deepcopy(PyMolConfig())

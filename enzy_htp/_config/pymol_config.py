"""Defines PyMOLConfig() which holds configuration settings for enzy_htp to interface with the 
PyMOL software package via the python interface. File also contains default_pymol_config() which
creates a default version of the PyMOLConfig() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-10-29
"""
from typing import Any, List
from copy import deepcopy


class PyMOLConfig:
    """Class that holds values for accessing the PyMOL API within enzy_htp. 
        
    Attributes:
        IO_EXTENSIONS: list() of str()'s with supported input/output file types.
    """
    #TODO(CJ): Info to store... which pymol version we are using


    IO_EXTENSIONS:List[str] = [".pdb", ".mol2", ".cif", ".pdb1"]
    """A list of the supported file extensions, including the leading dot. This is only a subset of pymol's supported filetypes."""


    def __getitem__(self, key: str) -> Any:
        """Getter that enables [] accession of PyMOLConfig() attributes."""
        if key.count("."):
            key1, key2 = key.split(".", 1)
            return getattr(self, key1)[key2]
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter that enables [] accession of PyMOLConfig() attributes with value validation."""
        if key.count("."):
            key1, key2 = key.split(".")
            AmberConfig.__dict__[key1][key2] = value
        else:
            setattr(self, key, value)

    
    def required_env_vars(self) -> List[str]:
        """There are no required environment variables for PyMOL at this time."""
        return list()


    def required_executables(self) -> List[str]:
        """There are no required executables for PyMOL at this time."""
        return list()

def default_pymol_config() -> PyMOLConfig:
    """Creates a deep-copied default version of the PyMOLConfig() class."""
    return deepcopy(PyMOLConfig())


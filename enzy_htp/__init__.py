"""
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-26
"""
# python imports
import sys
import logging
from typing import Tuple


# from .interface import Interface
from .structure import Structure, Residue
from .core import (
    EnvironmentManager,
    MissingEnvironmentElement,
    InvalidResidueCode,
    _LOGGER,
)
from .preparation import PDBLine, PDBPrepper, read_pdb_lines
from .chemical import ResidueType

from .molecular_mechanics import AmberInterface
from .quantum_mechanics import GaussianInterface, MultiwfnInterface

from .mutation import  MutationRestrictions, restriction_object, mutate_pdb


#CONFIG = None
#INTERFACE = None


# def init() -> Tuple[Config, Interface]:
#    """"""
#    global CONFIG
#    global INTERFACE
#    CONFIG = Config()
#    #INTERFACE = Interface(CONFIG)
#    return (CONFIG, INTERFACE)

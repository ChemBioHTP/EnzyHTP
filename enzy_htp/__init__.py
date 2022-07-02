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
	write_data
)
from .preparation import PDBLine, PDBPrepper, read_pdb_lines
from .chemical import ResidueType

from .molecular_mechanics import AmberInterface
from .quantum_mechanics import GaussianInterface, MultiwfnInterface

from .mutation import  MutationRestrictions, restriction_object, mutate_pdb

def welcome_msg():
    _LOGGER.info(f"""

#################################################### 
#      _____                 _   _ _____ ____      # 
#     | ____|_ __  _____   _| | | |_   _|  _ \     # 
#     |  _| | '_ \|_  / | | | |_| | | | | |_) |    #
#     | |___| | | |/ /| |_| |  _  | | | |  __/     # 
#     |_____|_| |_/___|\__, |_| |_| |_| |_|        # 
#                      |___/                       #
#                                                  # 
#################################################### 
""")


welcome_msg()

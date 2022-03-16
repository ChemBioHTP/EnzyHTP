# python imports
import sys
import logging

# enzy_htp imports
from .structure import Structure, Residue, PDBPrepper
from .core import EnvironmentManager, MissingEnvironmentElement, InvalidResidueCode, _LOGGER
from .chemical import ResidueType
from .config import CONFIG 


#def init_logging() -> None:
#    logging.basicConfig(
#        level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)], foj
#    )

#em.check_environment()

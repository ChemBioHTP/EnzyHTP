# python imports
import sys
import logging

# enzy_htp imports
from .structure import Structure, Residue, PDBPrepper
from .core import EnvironmentManager, MissingEnvironmentElement, InvalidResidueCode, _LOGGER
from .chemical import ResidueType

#def init_logging() -> None:
#    logging.basicConfig(
#        level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)], foj
#    )


em = EnvironmentManager(
    env_vars=["AMBERHOME"],
    executables=[
        "tleap",
        "ambpdb",
        "cpptraj",
        "mpirun",
        "$AMBERHOME/bin/sander.MPI",
        "$AMBERHOME/bin/pmemd.cuda",
        "$AMBERHOME/bin/MMPBSA.py.MPI",
    ],
)

em.check_environment()

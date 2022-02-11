""""""
import logging
from .exception import MissingEnvironmentElement
#tleap
#ambpdb
#cpptra
#mpirun
#$AMBERHOME
#$AMBERHOME/bin/sander.MPI
#$AMBERHOME/bin/pmemd.cuda
#$AMBERHOME/bin/MMPBSA.py.MPI


def __executable_in_path(exe_name: str) -> bool:
    pass


def __environment_variable_defined(var_name: str) -> bool:
    pass


def __compatible_os() -> str:
    pass


def check_compatability():
    # TODO need to setup the logging here
    logging.info("Beginning EnzyHTP compatability check...")
    error_ct = 0
    error_msg = []

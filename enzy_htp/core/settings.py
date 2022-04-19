"""Settings for the enzy_htp module. Responsibilities include versioning, path to data directory and determining operating system.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
2022-04-03
"""
from collections import namedtuple
import os
import sys

EnzyHTPVersion = namedtuple('EnzyHTPVersion', 'major minor patch')

_VERSION = _VERSION = EnzyHTPVersion(major = 0, minor = 1, patch = 0)

def get_version() -> EnzyHTPVersion:
    version_tuple = (_VERSION[0],_VERSION[1],_VERSION[2])
    return(version_tuple)

def version_to_str(version : EnzyHTPVersion) -> str:
    version_main = []
    for item in version:
        float_split = str(item).split('.')
        version_main.append(float_split[0])
    version_id = '.'.join(version_main)
    return(version_id)


def data_dir() -> str:
    #dirct = '../../enzy_htp/data/'
    #for i in os.listdir(dirct):
        #fulldirct = os.path.join(dirct,i)
        #full_dirct = fulldirct[5:]
        #return(full_dirct)
    dirct = 'enzy_htp/data/'
    return(dirct)

# returns path to a data/ directory in enzy_htp/data/

def is_compatible_os() -> bool:
    return sys.platform  in ['darwin', 'linux','linux2']
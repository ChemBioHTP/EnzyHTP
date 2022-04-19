"""Settings for the enzy_htp module. Responsibilities include versioning, path to data directory and determining operating system.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
2022-04-03
"""
from collections import namedtuple
import os.path
import sys

EnzyHTPVersion = namedtuple('EnzyHTPVersion', 'major minor patch')

_VERSION = EnzyHTPVersion(0,1,0)

def get_version() -> EnzyHTPVersion:  # should return version named tuple
    """Method that gets the current EnzyHTPVersion from the module."""
    return _VERSION

def version_to_str(version : EnzyHTPVersion) -> str:
    """Method that takes an EnzyHTPVersion and converts it to a str in the format "major.minor.patch"."""
    version = get_version()
    x = [ (int(a)) for a in version ]
    y = ".".join(str(elem) for elem in x)
    return y

def data_dir() -> str:
    """Returns a path to the data directory at enzy_htp/data/"""
    path = '/enzyhtp/data'
    return path

# returns path to a data/ directory in enzy_htp/data/

def is_compatible_os() -> bool:
    """Checks if the operating system is compatible. Should fail if operating system is not linux or mac-os"""
    if sys.platform != 'linux' or 'darwin':
        raise Exception('Compatible OS linux or macos only')

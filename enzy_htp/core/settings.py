"""Settings for the enzy_htp module. Responsibilities include versioning, path to data directory and determining operating system.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
2022-04-03
"""
from collections import namedtuple

EnzyHTPVersion = namedtuple('EnzyHTPVersion', 'major minor patch')

_VERSION = EnzyHTPVersion(major = 0, minor = 1,  patch = 0) #TODO(CJ): fill in the version here with major=0, minor=1 and patch = 0

def get_version() -> EnzyHTPVersion:  # should return version named tuple
    """Method that gets the current EnzyHTPVersion from the module."""
    return _VERSION

def version_to_str(version : EnzyHTPVersion) -> str:
    """Method that takes an EnzyHTPVersion and converts it to a str in the format "major.minor.patch".
	Should truncate floats:
	>>> version = EnzyHTPVersion(major=1.5, minor=3.2, patch=8.9)
	>>> version_to_str( version ) = "1.3.8"
	"""
    return str(int(version.major)) + '.' + str(int(version.minor)) + '.' + str(int(version.patch))

def data_dir() -> str:
    """Returns a path to the data directory at enzy_htp/data/"""
    return 'enzy_htp/data/'

# returns path to a data/ directory in enzy_htp/data/

def is_compatible_os() -> bool:
    """Checks if the operating system is compatible. Should fail if operating system is not linux or mac-os"""
	#HINT(CJ): Use sys.platform
    import sys
    system = sys.platform
    if system in ['linux', 'darwin']:
       return True
    else:
       return False

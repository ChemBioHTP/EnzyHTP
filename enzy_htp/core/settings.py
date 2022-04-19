"""Settings for the enzy_htp module. Responsibilities include versioning, path to data directory and determining operating system.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
2022-04-03
"""
from collections import namedtuple
import sys
EnzyHTPVersion = namedtuple('EnzyHTPVersion', 'major minor patch')

_VERSION = EnzyHTPVersion(0,1,0) #TODO(CJ): fill in the version here with major=0, minor=1 and patch = 0

def get_version() -> EnzyHTPVersion:  # should return version named tuple
    """Method that gets the current EnzyHTPVersion from the module."""
    return _VERSION

def version_to_str(version : EnzyHTPVersion) -> str:
    """Method that takes an EnzyHTPVersion and converts it to a str in the format "major.minor.patch".
	Should truncate floats:
	>>> version = EnzyHTPVersion(major=1.5, minor=3.2, patch=8.9)
	>>> version_to_str( version ) = "1.3.8"
	"""
    v_list = []
    s = map(lambda v: int(v), version)
    for x in  s:
        v_list.append(str(x))

    c = v_list[0]+'.'+v_list[1]+'.'+v_list[2]
    return c 

def data_dir() -> str:
    """Returns a path to the data directory at enzy_htp/data/"""
    return f'enzy_htp/data/'

# returns path to a data/ directory in enzy_htp/data/

def is_compatible_os() -> bool:
    """Checks if the operating system is compatible. Should fail if operating system is not linux or mac-os"""
	#HINT(CJ): Use sys.platform
    return sys.platform  in ['darwin', 'linux']

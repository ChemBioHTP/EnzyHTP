"""Information of the enzy_htp module. Responsibilities include versioning, path to data directory and determining operating system.
Should be import as the whole module if want to change constants. (e.g.: `from enzy_htp.core import enzyhtp_info`).

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-09-20
"""

from dataclasses import dataclass

@dataclass
class EnzyHTPVersion:
    """class that record version information of EnzyHTP
    We should follow versioning conventions in https://semver.org/"""
    major: int
    minor: int
    patch: int

    def __str__(self) -> str:
        """return the string representation of the version"""
        return f"{self.major}.{self.minor}.{self.patch}"

    def __lt__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

    def __le__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

    def __eq__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

    def __ne__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

    def __gt__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

    def __ge__(self, other) -> bool:
        raise Exception("this method is incomplele, complete this first!")

VERSION = EnzyHTPVersion(
        major=2,
        minor=0,
        patch=0,
    )

# TODO fill out the below functions

def get_version() -> EnzyHTPVersion:
    """get the current version of EnzyHTP."""
    return VERSION


# def data_dir() -> str:
#     pass


# # returns path to a data/ directory in enzy_htp/data/


# def is_compatible_os() -> bool:
#     pass
#     # checks if the operating system is compatbile

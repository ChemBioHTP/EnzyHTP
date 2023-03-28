"""Defines RosettaConfig() which holds configuration settings for enzy_htp to interface with the
Rosetta modelling suite. In addition to holding Rosetta. File also contains default_rosetta_config()
which creates a default version of the RosettaConfig() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""


from typing import List
from copy import deepcopy

class RosettaConfig:
    """ """

    pass
    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Rosetta."""
        return [
        ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Amber."""
        return [
        ]




def default_rosetta_config() -> RosettaConfig:
    """Creates a deep-copied version of the RosettaConfig() class."""
    return deepcopy(RosettaConfig())

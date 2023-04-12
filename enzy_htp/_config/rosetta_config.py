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

    ROSETTA3:str = "ROSETTA3"
    """ """

    ROSETTA_SCRIPTS:str = "rosetta_scripts.default.linuxgccrelease"#TODO(CJ): I think I need to change this
    """ """

    PARAMS_SCRIPT:str=f"ROSETTA3/source/scripts/python/public/molfile_to_params.py"


    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Rosetta."""
        return [
            self.ROSETTA_SCRIPTS,
            self.PARAMS_SCRIPT
        ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Amber."""
        return [
            self.ROSETTA3
        ]




def default_rosetta_config() -> RosettaConfig:
    """Creates a deep-copied version of the RosettaConfig() class."""
    return deepcopy(RosettaConfig())

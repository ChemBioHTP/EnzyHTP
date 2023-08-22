"""Defines RosettaConfig() which holds configuration settings for enzy_htp to interface with the
Rosetta modelling suite. In addition to holding Rosetta. File also contains default_rosetta_config()
which creates a default version of the RosettaConfig() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""

from typing import List, Any
from copy import deepcopy

from .base_config import BaseConfig

class RosettaConfig(BaseConfig):
    """Class that holds default values for running Rosetta with enzy_htp and 
    also keeps track of required environment variables and exectuables.

    Attributes:
        ROSETTA3 : str() of the environment variable needed for the quite.
        ROSETTA_SCRIPTS : str() of the RosettaScripts executable.
        PARAMS_SCRIPT: str() of the parameter generation file needed for RosettaLigand.
    """

    ROSETTA3: str = "ROSETTA3"
    """Path variable that points to the git installation of rosetta aka 'main'."""

    ROSETTA_SCRIPTS: str = "$ROSETTA3/source/bin/rosetta_scripts.linuxgccrelease"
    """The name of the RosettaScripts executable."""

    PARAMS_SCRIPT: str = f"$ROSETTA3/source/scripts/python/public/molfile_to_params.py"
    """Script used for paramterizing ligands for RosettaLigand protocol."""

    RELAX: str = f"$ROSETTA3/source/bin/relax.default.linuxgccrelease"
    """Executable used to relax a structure/pose."""

    SCORE: str = f"$ROSETTA3/source/bin/score_jd2.default.linuxgccrelease"
    """Executable used to score a specific structure/pose."""

    PY_2_7:str = "python2.7"
    """TODO(CJ)"""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Rosetta."""
        return [self.ROSETTA_SCRIPTS, self.PARAMS_SCRIPT, self.PY_2_7]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Rosetta."""
        return [self.ROSETTA3]

    def required_py_modules(self) -> List[str]:
        """ """
        return list()


def default_rosetta_config() -> RosettaConfig:
    """Creates a deep-copied version of the RosettaConfig() class."""
    return deepcopy(RosettaConfig())

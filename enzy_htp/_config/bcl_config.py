"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2023-04-02
"""


from copy import deepcopy


class BCLConfig:

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for the BCL."""
        return [
        ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for the BCL."""
        return []




def default_bcl_config() -> BCLConfig:
    """Creates a deep-copied default version of the BCLConfig() class."""
    return deepcopy(BCLConfig())

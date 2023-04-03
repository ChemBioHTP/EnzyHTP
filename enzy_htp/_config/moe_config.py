"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2023-04-02
"""


from copy import deepcopy


class MOEConfig:

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for the MOE."""
        return [
        ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for the MOE."""
        return []




def default_moe_config() -> MOEConfig:
    """Creates a deep-copied default version of the MOEConfig() class."""
    return deepcopy(MOEConfig())

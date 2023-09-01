""" """


from .base_config import BaseConfig

class RDKitConfig(BaseConfig):

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for xtb."""
        return list()

    def required_env_vars(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return list()


    def required_py_modules(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return ["rdkit"]

def default_xtb_config() -> XTBConfig:
    """Creates a deep-copied default version of the XTBConfig() class."""
    return deepcopy(XTBConfig())

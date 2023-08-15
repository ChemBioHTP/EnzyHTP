from copy import deepcopy
from .base_config import BaseConfig

class XTBConfig(BaseConfig):

    XTB_EXE:str = "xtb"

    def required_executables(self):
        return [self.XTB_EXE]


    def required_env_vars(self):
        return list()


    def required_py_modules(self):
        return list()



def default_xtb_config() -> XTBConfig:
    return deepcopy(XTBConfig())

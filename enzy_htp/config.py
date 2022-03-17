from .molecular_mechanics.amber_config import AmberConfig
from .molecular_mechanics.mmbpsa_config import MMPBSAConfig 
from .quantum_mechanics.gaussian_config import GaussianConfig
from .quantum_mechanics.multiwfn_config import MultiwfnConfig
from .core.env_manager import EnvironmentManager
from .core import _LOGGER
from copy import deepcopy

# TODO add documentation and logging for when things are changed

class Config:
    def __init__(self):
        self.ev_base =[]
        self.exe_base = ['tleap']
        self.env_vars = deepcopy(self.ev_base )
        self.executables = deepcopy(self.exe_base)
        self.n_cores = 24
        self.max_core = 2000
        self._amber_config = AmberConfig(self)
        self._mmpbsa_config = MMPBSAConfig(self)
        self._gaussian_config = GaussianConfig(self)
        self._multiwfn_config = MultiwfnConfig(self)
        self._em = None 
        self.update_paths()
        self.WORK_DIR = '.'

    def update_paths(self):
        self.env_vars = deepcopy(self.ev_base )
        self.executables = deepcopy(self.exe_base)

        self.executables.extend( self._amber_config.required_executables() )
        self.env_vars.extend( self._amber_config.required_env_vars() ) 
        
        self.executables.extend( self._mmpbsa_config.required_executables() )
        self.env_vars.extend( self._mmpbsa_config.required_env_vars() )
        
        self.executables.extend( self._gaussian_config.required_executables() )
        self.env_vars.extend( self._gaussian_config.required_env_vars() )
        
        self.executables.extend( self._multiwfn_config.required_executables() )
        self.env_vars.extend( self._multiwfn_config.required_env_vars() )
        self.executables = sorted(list(set(self.executables)))
        self.env_vars = sorted(list(set(self.env_vars)))
        
        self._em = EnvironmentManager(env_vars=self.env_vars,executables=self.executables)
        self._em.check_environment()

    def amber_config(self):
        return self._amber_config
    
    def mmbpsa_config(self):
        return self._mmbpsa_config
    
    def gaussian_config(self):
        return self._gaussian_config
    
    def gaussian_config(self):
        return self._gaussian_config

    def env_manager(self):
        return self._em

CONFIG = Config()

"""Defines a PyMOLInterface() class that serves as a bridge for enzy_htp to leverage the functionality 
provided by the python API of PyMOL. Uses the PyMOLConfig() class found in enzy_htp/_config/pymol_config.py.
Supported operations include charge calcuation. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-10-29
"""

from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em

from enzy_htp._config.pymol_config import PyMOLConfig, default_pymol_config

class PyMOLInterface:


    def __init__(self, config: PyMOLConfig = None ) -> None:
        self.config_ = config
        if not self.config_:
            self.config_ = default_amber_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            exectubles=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()

    def config(self) -> PyMOLConfig:
        """Getter for the PyMOLConfig() instance belonging to the class."""
        return self.config_



    def get_charge(self, fname : str, sele: str = 'all' ) -> int:
        

        #TODO(CJ): check that the correct file type is being used.
        pass

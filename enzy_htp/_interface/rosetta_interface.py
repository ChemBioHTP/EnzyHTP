"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em



from enzy_htp._config.rosetta_config import RosettaConfig, default_rosetta_config

class RosettaInterface:
    """ """



    def __init__(self, config: RosettaConfig = None ) -> None:
        """
        """
        self.config_ = config
        if not self.config_:
            self.config_ = default_amber_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            exectubles=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()


    def mutate(self):
        pass

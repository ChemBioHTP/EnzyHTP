from typing import List, Tuple
from pathlib import Path

from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.modeller_config import ModellerConfig, default_modeller_config




class ModellerInterface(BaseInterface):




    def __init__(self, parent, config: ModellerConfig = None ) -> None:
        super().__init__(parent, config, default_modeller_config)





    def fill_loops(self, ):
        pass


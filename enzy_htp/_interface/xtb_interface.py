



from .base_interface import BaseInterface

from enzy_htp._config.xtb_config import XTBConfig, default_xtb_config


class XTBInterface(BaseInterface):
    
    def __init__(self, parent, config) -> None:
        super().__init__(parent, config, default_xtb_config )




    def single_point(self, fname:str, charge:int=0, spin:int=1) -> float:

        """ """
        pass

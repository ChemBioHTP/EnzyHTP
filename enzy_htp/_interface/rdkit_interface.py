""" """



from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.rdkit_config import RDKitConfig, default_rdkit_config


class RDKitInterface(BaseInterface):
    pass


    def __init__(self, parent, config : RDKitConfig = None) -> None:
        """Simplicstic constructor that requires the parent interface as an argument and optionally takes an RDKitConfig instance.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_rdkit_config )


    def kekulize(self, molfile:str) -> str:
        pass

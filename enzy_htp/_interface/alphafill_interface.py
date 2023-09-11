"""
"""

from pathlib import Path

from .base_interface import BaseInterface



from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.alphafill_config import AlphaFillConfig, default_alphafill_config



class AlphaFillInterface(BaseInterface):
    

    def __init__(self, parent, config : AlphaFillConfig = None ) -> None:
        super().__init__(parent, config, default_alphafill_config )



    def fill_structure(self, molfile:str ) -> str:
        """
        """
        #TODO(CJ): add more options in here
        fs.check_file_exists( molfile )
        temp_path = Path(molfile)
        outfile = str(temp_path.parent / f"{temp_path.stem}_filled.cif")

        results = self.env_manager_.run_command(self.config_.ALPHAFILL_EXE,[
            "--config", self.config_.CONFIG_FILE,
            "process", molfile, outfile 
        ])


        return outfile

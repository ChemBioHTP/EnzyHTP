
from pathlib import Path

#TODO(CJ): add in itertions stuff
from .base_interface import BaseInterface

from enzy_htp._config.xtb_config import XTBConfig, default_xtb_config
from enzy_htp.core import file_system as fs

class XTBInterface(BaseInterface):
    
    def __init__(self, parent, config) -> None:
        super().__init__(parent, config, default_xtb_config )


    def single_point(self, fname:str, charge:int=0, spin:int=1) -> float:
        """ """
        results = self.env_manager_.run_command(self.config_.XTB_EXE,[
            f"--chrg", str(charge),
            "--sp",
            fname
        ])

        self._remove_temp_files(str(Path(fname).parent))

        for ll in results.stdout.splitlines():
            if ll.find('TOTAL ENERGY') != -1:
                return float(ll.split()[3])
        
        assert False


    def _remove_temp_files(self,work_dir:str) -> None:
        """
        """

        for fname in "wbo xtbrestart xtbtopo.mol".split():
            fs.safe_rm(f"{work_dir}/{fname}")

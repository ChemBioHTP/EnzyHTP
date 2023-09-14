"""Defines an AlphaFillInterface class that serves as a bridge for enzy_htp to utilize AlphaFill. Uses
the AlphaFillConfig class found in enzy_htp/_config/alphafill_config.py. Supported operations include:
    
    + filling structure with ligand transplants

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-14
"""
from pathlib import Path

from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.alphafill_config import AlphaFillConfig, default_alphafill_config


class AlphaFillInterface(BaseInterface):
    """ Class that provides a direct interface for enzy_htp to utilize AlphaFill. Supported operations
    include filling structures with ligand transplants. Users should use this class as the only way to 
    interact with this application.

    Attributes:
        config_ : The AlphaFillConfig() class which provides settings for both running AlphaFill and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """
    def __init__(self, parent, config : AlphaFillConfig = None ) -> None:
        """Simplistic constructor that optionally takes an AlphaFillConfig object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_alphafill_config )


    def fill_structure(self, molfile:str, outfile:str=None, work_dir:str=None, use_cache:bool=True ) -> str:
        """TODO(CJ): need to do this

        Args:
            molfile: The str() path to a .cif file to add ligand transplants to.
            outfile:
            work_dir:
            use_cache:
        
        Returns:
            The .cif file with transplanted ligands included.           
        """
        if fs.get_file_ext(molfile) != '.cif':
            _LOGGER.error(f"The supplied file {mofile} is not a .cif file! Exiting...")
            exit( 1 )

        #TODO(CJ): add more options in here
        fs.check_file_exists( molfile )
        temp_path = Path(molfile)
        outfile = str(temp_path.parent / f"{temp_path.stem}_filled.cif")
        if use_cache and fs.has_content(outfile):
            _LOGGER.info(f"The output file {outfile} exists and caching is enabled. Using this file.")
            return outfile
        else:
            fs.safe_rm( outfile )
                
        results = self.env_manager_.run_command(self.config_.ALPHAFILL_EXE,[
            "--config", self.config_.CONFIG_FILE,
            "process", molfile, outfile 
        ])

        return outfile

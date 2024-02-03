"""Defines an AlphaFillInterface class that serves as a bridge for enzy_htp to utilize AlphaFill. Uses
the AlphaFillConfig class found in enzy_htp/_config/alphafill_config.py. Supported operations include:
    
    + filling structure with ligand transplants
    + converting .pdb files into .mmcif files compatible with alphafill

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-14
"""
from typing import List, Tuple
from pathlib import Path

from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.alphafill_config import AlphaFillConfig, default_alphafill_config


class AlphaFillInterface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize AlphaFill. Supported operations
    include filling structures with ligand transplants. Users should use this class as the only way to 
    interact with this application.

    Attributes:
        config_ : The AlphaFillConfig() class which provides settings for both running AlphaFill and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: AlphaFillConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphaFillConfig object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_alphafill_config)


    def convert_file(self, infile:str, outfile:str=None) -> str:
        """Converts the supplied file to the specified output file. If no file is supplied, assumes mmCIF/.cif output. 
        When the input file is not a supported format, the function exits.

        Args:
            infile: Input file as a str().
            outfile: File to convert the input to. Optional. When not specified creates the input file with a .cif extension. 

        Returns:
            The converted mmCIF file.
        """
        if outfile is None:
            outfile = str(Path(infile).with_suffix('.cif'))            

        _MAPPER:Dict[Tuple[str,str],str] = { (".pdb", ".cif"):"1", (".cif", ".pdb"):"2", (".cif", ".cif"):"8" }
        start_ext:str=fs.get_file_ext(infile)
        end_ext:str=fs.get_file_ext(outfile)

        conversion_code:str=_MAPPER.get((start_ext, end_ext), None)

        if conversion_code is None:
            _LOGGER.error(f"The conversion from {start_ext} to {end_ext} is not supported! Exiting...")
            exit( 1 )

        self.env_manager_.run_command(self.config_.MAXIT , [
            "-input", infile, "-output", outfile, "-o", conversion_code
        ] ) 

        fs.safe_rm( "./maxit.log" )

        return outfile

    def fill_structure(self, molfile: str, outfile: str = None, work_dir: str = None, use_cache: bool = True) -> Tuple[str, str]:
        """Using MSA-derived constraints, place ligands into the supplied structure file. File can be either .cif or .pdb format.
        Note that if it is .cif format it MUST be mmCIF, not just .cif

        Args:
            molfile: The str() path to a .cif file to add ligand transplants to.
            outfile: Where the filled structure will be saved to. Optional.
            work_dir: Where temporary files will be saved to. Optional.
            use_cache: Should we use existing files when available? Optional.
        
        Returns:
            A Tuple[str, str] with the format ( .cif file with transplants, .json file from alphafill). 
        """
        if work_dir is None:
            work_dir = self.parent().config()['system.SCRATCH_DIR']

        fs.safe_mkdir( work_dir )

        if fs.get_file_ext(molfile) == ".pdb":
            _LOGGER.info(f"The supplied file {molfile} is .pdb format. Converting to mmCIF...")
            molfile = self.convert_file(molfile)

        if fs.get_file_ext(molfile) != '.cif':
            _LOGGER.error(f"The supplied file {mofile} is not a .cif file! Exiting...")
            exit(1)

        #TODO(CJ): add more options in here
        fs.check_file_exists(molfile)
        temp_path = Path(molfile)
        outfile = str(temp_path.parent / f"{temp_path.stem}_filled.cif")
        json_outfile = str(temp_path.parent / f"{temp_path.stem}_filled.json")
        
        if use_cache and fs.has_content(outfile) and fs.has_content( json_outfile ):
            _LOGGER.info(f"The output files {outfile} and {json_outfile} exist and caching is enabled. Using these files as is.")
            return (outfile, json_outfile)
        else:
            fs.safe_rm(outfile)

        results = self.env_manager_.run_command(self.config_.ALPHAFILL_EXE,
                                                ["--config", self.config_.CONFIG_FILE, "process", molfile, outfile])
        fs.safe_rm(temp_path)
        return (outfile, json_outfile)

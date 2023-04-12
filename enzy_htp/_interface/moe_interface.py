"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-04-01
"""

from pathlib import Path
from typing import List, Dict

from ..core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em

from enzy_htp._config.moe_config import MOEConfig, default_moe_config

#TODO(CJ): need to add tests for this section

class MOEInterface:
    
    def __init__(self, config: MOEConfig = None) -> None:
        """Simplistic constructor that optionally takes an MOEConfig object as its only argument.
        Also checks if the current environment is compatible with the MOEInterface().
        """
        self.config_ = config
        if not self.config_:
            self.config_ = default_moe_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            executables=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()


    def valid_license(self) -> bool:
        """MOE operates using a license based model. This function checks if a valid license is available in the system."""
        pass


    def _run_moebatch(self, input_svl:str) -> None:
        """Helper function that runs the supplied scientific vector language (SVL) file as input with moebatch.
        Method checks that the supplied svl file exists and has the correct file extension.

        Args:
            input_svl: A str() with the path of the .svl file.

        Returns:
            Nothing
        """
        
        fs.check_file_exists( input_svl )

        if Path(input_svl).suffix != '.svl':
            _LOGGER.error(f"The supplied svl file '{input_svl}' does not have the extension .svl. Exiting...")
            exit( 1 )


        self.env_manager_.run_command(
            self.config_.MOE_BATCH,
            ["-s", "input_svl"]
        )


    def protonate(self, molfile:str, outfile:None=str, pH:float=7.0) -> str:
        """Protonates the structure in a supplied file at the specified pH. If nothing is supplied for the
        outfile, then it is set as <input_file>_protonated.<ext> for a given file of name <input_file>.<ext>.
        Checks that the supplied molfile exists and that it has one of the supported file types of .mol2,
        .pdb, or .sdf. Also checks that the supplied pH is on the range [0.0, 14.0]. If these criteria are
        not met, the method exits.

        Args:
            mofile: File containing the input structure. Acceptable formats include .mol2, .pdb and .sdf.
            outfile: File name for the protonated molecule. Optional argument.
            pH: A float() with the pH to protonate the ligand at. Default is 7.0.

        Returns:
            The path of the protonated structure as a str().            
        """
        
        fs.check_file_exists( molfile )
        
        FORMAT_MAPPER:Dict[str,str] = {
            '.mol2':'TriposMOL2',
            '.pdb':'PDB',
            '.sdf':'SDF',
        }
        
        ftype:str = FORMAT_MAPPER.get(Path(molfile).suffix, None)

        if not outfile:
            tpath = Path( molfile )
            outfile = tpath.parent / f"{tpath.stem}_protonated{tpath.suffix}"

        if not ftype:
            _LOGGER.error(f"Supported file types for MOEInterface.protonate() only include .mol2, .pdb, and .sdf. Exiting...")
            exit( 1 )

        #TODO(CJ): need to figure out where I modify the pH
        #TODO(CJ): check that the pH value is valid/accurate
        contents:List[str]=[
            f"ReadTriposMOL2 '{molfile}';", 
             "atoms = Atoms[];",
             "Protonate3D [atoms,atoms,atoms,[],[],[]];",
            f"WriteTriposMOL2 '{outfile}';",
             "Close [];"
        ]
       
        temp_file:str = Path(molfile).parent / "temp.svl"
        
        fs.write_lines( temp_file, contents )

        self._run_moebatch( temp_file )
        
        fs.safe_rm( temp_file )

        return outfile


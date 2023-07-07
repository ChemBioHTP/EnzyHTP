"""Defines a BCLInterface class that serves as a bridge for enzy_htp to utilize the biochemical library (BCL)
Uses the BCLConfig class found in enzy_htp/_config/_config.py. Supported operations include conformer generation.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-04-01
"""
import pandas as pd
from typing import List
from pathlib import Path

from ..core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em

from .base_interface import BaseInterface

from enzy_htp._config.bcl_config import BCLConfig, default_bcl_config


class BCLInterface(BaseInterface):
    """Class that provides a direct interfrace for enzy_htp to utilize the BCL. Supported opterations
    include conformer generation. Users should use this class as the only way to interact with any application
    in BCL.


    Attributes:
        config_	: The BCLConfig() class which provides settings for both running BCL and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensure all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: BCLConfig = None) -> None:
        """Simplistic constructor that optionally takes a BCLConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_bcl_config)

    def _run_bcl_cmd(self, protocol: str, flags: List[str]) -> None:
        """Underlying method to use the BCL. BCL has a format/calling convention of:
            $ bcl.exe <protocol> <flags>
        
        This function is NOT meant to be used directly but should be used elsewhere in the BCLInterface class.

        Args:
            protocol: The name of the BCL protocol to use as a str(). 
            flags: A list() of str() flags to pass to bcl.exe.

        Returns:
            Nothing.

        """
        flags = [protocol] + flags
        self.env_manager_.run_command(self.config_.BCL_EXE, flags)

    def generate_conformers(self, molfile: str, n_conformers: int = 100, outfile: str = None) -> str:
        """Method that creates conformers using the 'molecule:ConformerGenerator' method from the BCL.
        Note that the supplied molfile must exist but outfile does not. All conformers are saved to one
        file. When the outfile is not supplied, the output file is saved to the same directory as the 
        input file with the suffix '_conformers' added to the name.

        Args:
            molfile: str() with the name of the input file containing the small molecule of interest. 
            n_conformers: The int() number of conformers to make.
            outfile: str() with the path to save the conformers to. Optional argument.


        Returns:
           A str() naming the file with the generated conformers. 
        """

        #TODO(CJ): add the ability to add more flags and that molfile is in the correct format

        fs.check_file_exists(molfile)

        if not outfile:
            fpath = Path(molfile)
            outfile = fpath.parent / f"{fpath.stem}_conformers.sdf"

        flags: List[str] = [
            "-conformation_comparer",
            "SymmetryRMSD",
            "0.25",
            "-ensemble_filenames",
            molfile,
            "-top_models",
            f"{int(n_conformers)}",
            "-conformers_single_file",
            f"{outfile}",
        ]

        self._run_bcl_cmd("molecule:ConformerGenerator", flags)

        return outfile

    def formal_charge(self, molfile: str) -> int:
        """Find the formal charge of the supplied molfile using molecule:Properties routine. Supplied
        file MUST be in format .sdf otherwise the script will exit.

        Args:
            molfile: A str() with the path of an .sdf file.

        Returns:
            Formal charge of the molecule as an int().
        """

        fs.check_file_exists(molfile)

        if not Path(molfile).suffix == '.sdf':
            _LOGGER.error(f"Function expects .sdf file format. Exiting...")
            exit(1)

        temp_file: str = "__temp.csv"

        flags: List[str] = ["-input_filenames", str(molfile), "-output_table", temp_file, "-tabulate", "TotalFormalCharge"]

        self._run_bcl_cmd("molecule:Properties", flags)

        df: pd.DataFrame = pd.read_csv(temp_file)

        fs.safe_rm(temp_file)

        return df.iloc[0].TotalFormalCharge

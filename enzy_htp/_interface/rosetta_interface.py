"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
from typing import List
from pathlib import Path

import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em



from enzy_htp._config.rosetta_config import RosettaConfig, default_rosetta_config

class RosettaInterface:
    """Class that interfaces with the Rosetta software package.


    """



    def __init__(self, config: RosettaConfig = None ) -> None:
        """Trivial constructor for the RosettaInterface class. Can take no arguments or a RosettaConfig object.

        Args:
            config: Optional. A RosettaConfig 
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


    def run_rosetta_scripts(self, pdb:str, script_file:str, opts:List[str]) -> None:
        """Method that runs the rosettascripts executable, allowing for  
        """

        options:List[str] = [
            "-si",
            str(pdb),
            "-parser:protocol",
            str(script_file)# -overwrite -mute all -out:prefix work_dir/
        ]
        
        options.extend( opts )

        self.env_manager_.run(
            self.config_.ROSETTA_SCRIPTS,
            opts
        )


    def parse_score_file(self, fname : str) -> pd.DataFrame:
        """Method that parses a score file into a Pandas Dataframe. Only encodes lines that begin with SCORE.

        Args:
            fname: Path to the score file. Will error if does not exist.

        Returns:
            A pandas dataframe containing the data in the supplied score file.
        """
        
        if not Path(fname).exists():
            _LOGGER.error(f"The suppliied file '{fname}' does not exist. Exiting...")
            exit( 1 )
        
        lines:List[str] = fs.lines_from_file( fname )
        
        lines = list(filter(lambda ll: ll.startswith('SCORE:'), lines ))
        
        data = list(map(lambda ll: ll.split()[1:], lines ))

        df = pd.DataFrame(data=data[1:], columns=data[0])
       
        column_names = list(df.columns)
        for cn in column_names:
            if cn == 'description':
                continue
            df[cn] = df[cn].astype('float')

        return df

#/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease  #TODO(CJ): this rosetta_scripts exe does not exist in this form

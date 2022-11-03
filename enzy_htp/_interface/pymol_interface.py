"""Defines a PyMOLInterface() class that serves as a bridge for enzy_htp to leverage the functionality 
provided by the python API of PyMOL. Uses the PyMOLConfig() class found in enzy_htp/_config/pymol_config.py.
Supported operations include charge calcuation. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-10-29
"""
from pathlib import Path

from typing import List, Dict, Any

from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em
from enzy_htp.core import _LOGGER

from enzy_htp._config.pymol_config import PyMOLConfig, default_pymol_config

class PyMOLInterface:
    """
    """

    def __init__(self, config: PyMOLConfig = None ) -> None:
        from pymol import cmd
        self.cmd = cmd
        self.config_ = config
        if not self.config_:
            self.config_ = default_amber_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            exectubles=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()

    def config(self) -> PyMOLConfig:
        """Getter for the PyMOLConfig() instance belonging to the class."""
        return self.config_


    def supported_file_type(self, fname : str ) -> bool:
        """Convenience function that checks if a listed filetype is supported by 
        the PyMOLInterface. Currently supported filetypes include:
            + .pdb
            + .mol2
            + .cif

        Args:
            fname: The str() path to the file in question.

        Returns:
            Whether the file type is supported.
        """
        extension:str = Path(fname).suffix
        return extension in self.config().IO_EXTENSIONS

    def get_charge(self, fname : str, sele: str = '(all)' ) -> int:
        """Method that gets the formal charge for the specified sele in the 
        specified file. File must be a supported file type as listed in 
        PyMOLConfig. Checks if file exists and is supported type. If either 
        are not true then an error is logged and the program exits. DOES
        NOT check if the sele is valid. 

        Args:
            fname: The str() path to the file in question.
            sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.

        Returns:
            The formal charge of the sele as an int().
        """
        if not Path(fname).exists():
            _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
            exit( 1 )

        if not self.supported_file_type( fname ):
            _LOGGER.error(f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit( 1 )

        _eh_local : Dict[str, Any] = { 'fc': 0 }


        self.cmd.delete('all')
        self.cmd.load( fname )
        self.cmd.iterate( sele, 'fc += formal_charge', space=_eh_local )
        self.cmd.delete('all')


        return _eh_local['fc']
        

    def get_sequence(self, fname : str, sele: str = '(all)' ) -> str:
        """Method that gets the sequence for the specified sele in the 
        specified file. File must be a supported file type as listed in 
        PyMOLConfig. Checks if file exists and is supported type. If either 
        are not true then an error is logged and the program exits. DOES
        NOT check if the sele is valid. Note that non-canonical residues will
        be represented by a '?' in the sequence.

        Args:
            fname: The str() path to the file in question.
            sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.

        Returns:
            A str() with the sequence of the sele as it would appear in a .fasta file..
        """
        if not Path(fname).exists():
            _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
            exit( 1 )

        if not self.supported_file_type( fname ):
            _LOGGER.error(f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit( 1 )
       
        self.cmd.delete('all')
        self.cmd.load( fname )
        lines:List[str] = self.cmd.get_fastastr( sele ).splitlines()
        self.cmd.delete('all')

        lines = list(filter( lambda ll: ll[0] != '>', lines))
        return ''.join(lines)

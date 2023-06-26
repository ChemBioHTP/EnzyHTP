"""Defines a PyMOLInterface() class that serves as a bridge for enzy_htp to leverage the functionality 
provided by the python API of PyMOL. Uses the PyMOLConfig() class found in enzy_htp/_config/pymol_config.py.
Supported operations include charge calcuation. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-10-29
"""
from pathlib import Path

from typing import List, Dict, Any, Union, Tuple

import pandas as pd

from enzy_htp.core import _LOGGER, check_var_type
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em

from enzy_htp._config.pymol_config import PyMOLConfig, default_pymol_config

from .base_interface import BaseInterface

#TODO(CJ): add something to remove "PyMOL not running. Entering library mode (experimental" message on pymol running

class PyMOLInterface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize the PyMOL python package. Supported operations 
    include (TODO(CJ))


    Attributes:
        cmd : The pymol module that enzy_htp as has ccess to.
        config_ : The PyMOLConfig() class which provides settings for both running PyMOL and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """
    def __init__(self, parent, config: PyMOLConfig = None) -> None:
        """Simplistic constructor that optionally takes an PyMOLConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_pymol_config)
        #TODO(CJ): set it up so I can import pymol OR pymol2
        import pymol 
        pymol.invocation.parse_args(['pymol', '-q'])

        from pymol import cmd
        self.cmd = cmd
        self.available_cmds_:List[str] = dir( self.cmd )

    def convert(self, file_1:str, file_2:str=None, new_ext:str=None, split_states:bool=False) -> Union[str, List[str]]:
        """Method that converts a supplied file to a different format. Either a new filename or new file 
        extension can be supplied. If neither or both are supplied, then the function will exist. Note
        that the function does not check for valid file types and will catch any errors that are thrown
        if an invalid file_1 or output file combination is supplied. Returns the new outfile.

        Args:
            file_1 : The name of the original file as a str().
            file_2 : The name of the output file as a str(). Optional.
            new_ext : The new extension to use. Optional.

        Returns:
            The name of the new file as a str().
        """
        #TODO(CJ): update split states. update return policy
        
        fs.check_file_exists( file_1 )

        if (not file_2 and not new_ext) :
            _LOGGER.error("Either a new file name or target extension must be supplied. Both are empty. Exiting...")
            exit( 1 )

        if file_2 and new_ext:
            _LOGGER.error("Either a new file name or target extension must be supplied. Both were supplied. Exiting...")
            exit( 1 )

        if not file_2:
            fpath = Path( file_1 )
            file_2 = fpath.with_suffix( new_ext )
        
        result:str = str(file_2)

        try:
            self.cmd.delete('all')
            self.cmd.load(file_1)
            
            if split_states:
                self.cmd.split_states('all')                
                result = Path(result)
                parent_dir = result.parent 
                result_stem = result.stem
                result = []
                
                for oidx, oname in enumerate(self.cmd.get_object_list()):
                    
                    if self.cmd.count_states(oname) != 1:
                        self.cmd.delete(oname)
                        continue
                    
                    fname = str( parent_dir / f"{result_stem}_{oidx}{new_ext}" )
                    self.cmd.save( fname, oname )
                    result.append( fname )
            else:
                self.cmd.save(result)
                self.cmd.delete('all')
        
        except Exception as exe:
            _LOGGER.error(f"Could not convert '{file_1}' to '{file_2}'. Encountered error: '{exe}'. Exiting...")
            exit( 1 )

        return result 




    def supported_file_type(self, fname: str) -> bool:
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
        #TODO(CJ): need to check which file formats are actually supported
        extension: str = Path(fname).suffix
        return extension in self.config().IO_EXTENSIONS

    def get_charge(self, fname: str, sele: str = '(all)') -> int:
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
            exit(1)

        if not self.supported_file_type(fname):
            _LOGGER.error(
                f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit(1)

        _eh_local: Dict[str, Any] = {'fc': []}

        self.cmd.delete('all')
        self.cmd.load(fname)
        self.cmd.iterate(sele, 'fc.append(formal_charge)', space=_eh_local)
        self.cmd.delete('all')

        return sum(_eh_local['fc'])

    def get_sequence(self, fname: str, sele: str = '(all)') -> str:
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
            exit(1)

        if not self.supported_file_type(fname):
            _LOGGER.error(
                f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit(1)

        self.cmd.delete('all')
        self.cmd.load(fname)
        lines: List[str] = self.cmd.get_fastastr(sele).splitlines()
        self.cmd.delete('all')

        lines = list(filter(lambda ll: ll[0] != '>', lines))
        return ''.join(lines)

    def de_protonate(self, molfile : str, new_file:str=None) -> str:
        """Method that deprotonates the structure supplied by molife, saving it to a new file. 
        If no new_file is supplied, original file <name>.<ext> will be saved to <name>_deprotonated.<ext>

        Args:
            molfile: The name of the input file.
            new_file: Where to save the deprotonated structure. Optional.
    
        Returns:
            The name of the deprotonated structure.
        """

        fs.check_file_exists( molfile )
        
        if not new_file:
            fpath = Path(molfile )
            new_file = fpath.parent / f"{fpath.stem}_deprotonated{fpath.suffix}"

        self.cmd.delete( 'all' )
        self.cmd.load( molfile )
        self.cmd.remove( 'hydrogens' )
        self.cmd.save( new_file )

        self.cmd.delete('all')

        return new_file

    def rename_atoms(self, template: str, in_file: str) -> str:
        """ """
        _eh_local: Dict[str, Any] = {'template': [], 'in_stru':[]}

        self.cmd.delete('all')
        self.cmd.load(template)
        self.cmd.iterate('all', 'template.append((name, elem, ID))', space=_eh_local)
        self.cmd.delete('all')

        self.cmd.load(in_file)
        self.cmd.iterate('all', 'in_stru.append((name, elem, ID))', space=_eh_local)

        if len(_eh_local['template']) != len(_eh_local['in_stru']):
            _LOGGER.error(f"Atom number mismatch in template vs in_file structures. ({len(_eh_local['template'])} vs {len(_eh_local['in_stru'])}). Exiting...")
            exit( 1 )
        
        for tidx, (tp1, tp2) in enumerate(zip(
                sorted(_eh_local['template'], key=lambda tp: tp[-1]),
                sorted( _eh_local['in_stru'], key=lambda tp: tp[-1])
                 )):
            if tp1[1] != tp2[1]:
                _LOGGER.error(f"Atom type mismatch at atom number {tidx+1} template vs in_file structures. template: {tp1[1]}  in_stru: {tp2[1]}. Exiting...")
                exit( 1 )

        for (name, _, ID) in _eh_local['template']:
            self.cmd.alter( f"ID {ID}", f'name="{name}"' )

        self.cmd.save( in_file )
        self.cmd.delete('all')

        return in_file


    def collect(self, molfile:str, variables:List[str], sele:str='all', state:int=-1) -> pd.DataFrame:
        """Function that serves as a wrapper to the pymol cmd API's iterate_state() method. Checks if the supplied file exists prior to collecting data.
        Data is returned in a pandas DataFrame() where each row is an atom and each column is a variable.

        Args:
            molfile: Name of the structure to sample. 
            variables: A list() of str()'s with variable names to pull from in the structure.
            sele: The pymol selection string to use. Defaults to 'all'
            state: The state to use. Defaults to -1 (current state).

        Returns:
            A pandas DataFrame() where each row is an atom and each column is a variable.

        """


        check_var_type(variables, list)
        for vv in variables:
            check_var_type( vv, str )

        if molfile  != 'memory':
            fs.check_file_exists(molfile)    
            self.cmd.delete('all')

        _eh_local: Dict[str, Any] = {'data':[]}

        if molfile != 'memory':
            self.cmd.delete('all')
            self.cmd.load(molfile)

        iter_stmt:str=f"data.append(({', '.join(variables)}))"
        self.cmd.iterate_state(state,sele,iter_stmt,space=_eh_local)

        result = pd.DataFrame(
            data=_eh_local['data'],
            columns=variables
        )

        if molfile != 'memory':
            self.cmd.delete('all')
        
        return result

    def execute(self, args:List[Tuple]) -> List[Any]:
        """Executes a series of commands through the PyMOL/PyMOL2 python module in use. Takes input as a list of Tuple()'s
        where the first item in each tuple is a string specifying the function to use and the rest of the items are the
        arguments for that function.
    
        Args:
            args: A list() of tuple()'s of items where the first string names the pymol API command to use and the rest are arguments for that function.

        Returns:
            Each result for the respective command output.
        """
       
        result:List[Any] = list() 

        for cmd_set in args:
            if len(cmd_set) < 2:
                _LOGGER.error(f"The supplied argument {cmd_set} is not long enough. Musst have at least two items. Exiting...")
                exit( 1 )

            cmd_name = cmd_set[0]
            cmd_args = list(cmd_set[1:])
            
            if cmd_name not in self.available_cmds_:
                _LOGGER.error(f"The command '{cmd_name}' is not supported in this version of pymol. Exiting...")
                exit( 1 )
           

            cmd_str:str=f"{cmd_name}({','.join(map(str, cmd_args))})"
            try:
                result.append(
                    self.cmd.__dict__[cmd_name](*cmd_args)
                )
            except:
                _LOGGER.error(f"PyMOL function call '{cmd_str}' resuled in an error. Exiting...")
                exit( 1 )





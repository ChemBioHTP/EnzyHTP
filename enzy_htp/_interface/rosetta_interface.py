"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
import shutil
from pathlib import Path
from typing import List, Tuple
from collections import namedtuple

import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em



from enzy_htp._config.rosetta_config import RosettaConfig, default_rosetta_config


#TODO(CJ): make sure the list() of str()'s are only strings

RosettaResult=namedtuple(
    'RosettaResult',
    'score_sc structures log_lines'
)
RosettaResult.__doc__='TODO(CJ)'


from .base_interface import BaseInterface

class RosettaInterface(BaseInterface):
    """Class that interfaces with the Rosetta software package.

    Attributes:
        

    """

    def __init__(self, parent, config: RosettaConfig = None) -> None:
        """Simplistic constructor that optionally takes an RosettaConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_rosetta_config)



    def run_rosetta_scripts(self, opts:List[str], logfile:str=None) -> None:
        """Method that runs the rosettascripts executabl along with the supplied options. Optionally outputs
        the stdout to a specified logfile. Note that no sanitation is performed prior to running.

        Args:
            opts: a list() of str() to be run by the RosettaScripts executable.
            logfile: The file to output the stdout log to. Optional.

        Returns:
            Nothing,
        """

        if logfile:
            opts.extend(
                [
                    ">",str(logfile)
                ]
            )

        self.env_manager_.run_command(
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


    def parameterize_ligand(self, molfile:str, res_name:str, outdir:str=None, conformers:str=None) -> Tuple[str, str]:
        """Parameterizes the input ligand for use in the RosettaLigand protocol. Takes an input file with the ligand,
        as well as the name of the residue in PDB format (3 capitalized letters) as well as optionally where the output
        directory where the .params and .pdb files should be saved. The underlying script only supports .mol, .mol2, 
        and .sdf formats. Can also add conformers to end of .params file when conformers file is supplied. Function 
        exits on invalid inputs. 

        Args:
            molfile: The name of the input file as a str().
            res_name: The all-capitalized, three letter string of the ligand in PDB format.
            outfile: Where the .params file should be saved. Optional.
            conformers: The conformers file for the given ligand. Optional.

        Returns:
            A tuple() with the layout of (.params file, .pdb file)            

        """
        #TODO(CJ): add in the ability to add a conformers file
        fs.check_file_exists( molfile )
    
        ALLOWED_FORMATS:List[str] = ".sdf .mol2 .mol".split()
        
        if Path(molfile).suffix not in ALLOWED_FORMATS:
            _LOGGER.error(f"The supplied file '{molfile}' is of an unsupported file type. Supported types are {', '.join(ALLOWED_FORMATS)}. Exiting...")
            exit( 1 )
        #TODO(CJ): need to fix this. codes like 152 are valid  => only numbers
        #if not res_name.isupper() or not len(res_name) == 3 or not res_name.isalpha():
        #    _LOGGER.error(f"The supplied residue name '{res_name}' is invalid. It must be alphanumeric, capitalized, and have three characters. Exiting...")
        #    exit( 1 )
        
        flags:List[str] = [
            f"{molfile}",
            f"--name={res_name}",
            "--clobber"
        ]

        self.env_manager_.run_command(
            self.config_.PARAMS_SCRIPT,
            flags
        )

        params_file:str = f"./{res_name}.params"
        pdb_file:str = f"./{res_name}_0001.pdb"

        if outdir:
            fs.safe_mkdir( outdir )
            
            outdir = Path(outdir)
            params_start = params_file
            pdb_start = pdb_file
            
            params_file = str( outdir / Path(params_start).name)
            pdb_file = str(outdir / Path(pdb_start).name)

            shutil.move(params_start, params_file)
            shutil.move(pdb_start, pdb_file)

        fs.check_file_exists( params_file )
        fs.check_file_exists( pdb_file )


        if conformers:
            self.add_conformers( params_file, conformers )

        return (params_file, pdb_file)


    def add_conformers(self, param_file:str, conformers_file:str) -> None:
        """Adds conformers to the end of the .params file for a given ligand. Checks that files exist and
        exits if not.

        Args:
            param_file: The .params file as a str().
            conformers_file: The path to the conformers file. MUST BE IN .pdb FORMAT. 

        Returns:
            Nothing.
        """
        
        fs.check_file_exists( param_file )
        fs.check_file_exists( conformers_file )

        if not Path(conformers_file).suffix == ".pdb":
            _LOGGER.error(f"The supplied file '{conformers_file}' is not in the .pdb format. Exiting...")
            exit( 1 )

        content:List[str] = fs.lines_from_file( param_file )
        content.append(f"PDB_ROTAMERS {Path(conformers_file).absolute()}") 
        fs.write_lines(param_file, content)


    def relax(self, 
        infile:str, 
        nstruct:int, 
        ignore_zero_occupancy:bool=True,
        full_atom:bool=True,
        detect_disulf:bool=True,
        linmem_ig:int=10,
        constrain_relax_to_start_coords:bool=True,
        coord_constrain_sidechains:bool=True, 
        ramp_constraints:bool=True,
        prefix:str=None,
        overwrite:bool=True,
        extra_flags:List[str]=None
        #work_dir:str=None
        ) -> str:
        """ """
        #TODO(CJ): need to add documentation; should also be able to take a Structure as input.
        fs.check_file_exists( infile )

            #/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/relax.default.linuxgccrelease
            #-out:prefix $prefix
            #-out:file:scorefile ${prefix}.sc &
        flags:List[str]=[
            f"-in:file:s '{infile}'",
            f"-nstruct {nstruct}",
            f"-linmem_ig {linmem_ig}",
        ]
        
        
        flags.append(f"-ignore_zero_occupancy {'true' if ignore_zero_occupancy else 'false'}")
        flags.append(f"-relax:constrain_relax_to_start_coords {'true' if constrain_relax_to_start_coords else 'false'}")
        flags.append(f"-coord_constrain_sidechains {'true' if coord_constrain_sidechains else 'false'}")
        flags.append(f"-ramp_constraints {'true' if ramp_constraints else 'false'}")


        if full_atom:
            flags.append("-in:file:fullatom")

        if detect_disulf:
            flags.append("-in:detect_disulf")

        if prefix:
            flags.append(f"-out:prefix '{prefix}'")

        if overwrite:
            flags.append("-overwrite")

        if extra_flags:
            flags.extend(extra_flags)

        self.env_manager_.run_command(
            self.config_.RELAX,
            flags
        )
        #TODO(CJ): what do I want to return from here?


    def score(self,
            infile:str,
            ignore_zero_occupancy:bool=True,
            overwrite:bool=True,
            extra_flags:List[str]=None
        ) -> float:
        """Provides the total score in Rosetta Energy Units (REU) for a given structure. Uses default flags but can have behavior modified
        via supplied extra_flags. Returns the total score in REU.

        Arguments:
            infile: The file to 
            ignore_zero_occupancy: 
            overwrite:
            extra_flags:

        Returns:
            Score of structure in file in REU.

        """
        fs.check_file_exists( infile )
        
        flags:List[str]=[
            f"-in:file:s '{infile}'",
             "-ignore_unrecognized_res",
        ]
        
        flags.append(f"-ignore_zero_occupancy {'true' if ignore_zero_occupancy else 'false'}")

        if overwrite:
            flags.append("-overwrite")

        if extra_flags:
            flags.extend(extra_flags)

        fs.safe_rm( "./score.sc" )

        self.env_manager_.run_command(
            self.config_.SCORE,
            flags
        )


        df:pd.DataFrame=self.parse_score_file("./score.sc")
    
        #TODO(CJ): figure this out/make it better
        assert len(df) == 1

        return df.iloc[0].total_score


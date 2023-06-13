"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation,
relaxation (minimization), scoring, ligand parameterization, and the ability to use RosettaScripts.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
import shutil
from pathlib import Path
from collections import namedtuple
from typing import List, Tuple, Dict

import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em



from enzy_htp._config.rosetta_config import RosettaConfig, default_rosetta_config

from .base_interface import BaseInterface

class RosettaInterface(BaseInterface):
    """Class that provides an interface for enzy_htp to utilize the Rosetta molecular modelling software suite. Supported 
    operations include relaxation (minimzation), scoring, ligand parameterization, and the ability to use RosettaScripts.

    Attributes:
        config_	: The RosettaConfig() class which provides settings for both running Rosetta and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: RosettaConfig = None) -> None:
        """Simplistic constructor that optionally takes an RosettaConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_rosetta_config)

    def _delete_score_file(self, working_dir:str='./') -> None:
        """Helper method that deletes the score.sc file in the specified directory.

        Args:
            working_dir: The directory to look for the score.sc file in. Defaults to './'

        Returns:
            Nothing.
        """
        score_file = working_dir + "/score.sc"
        fs.safe_rm( score_file )


    def _delete_crash_log(self, working_dir:str='./') -> None:
        """Helper method that deletes the ROSETTA_CRASH.log file in the specified directory.

        Args:
            working_dir: The directory to look for the ROSETTA_CRASH.log file in. Defaults to './'

        Returns:
            Nothing.
        """
        log_name = working_dir + "/ROSETTA_CRASH.log"
        fs.safe_rm( log_name )

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

        if logfile:
            _LOGGER.info(f"Saved RosettaScripts log to '{Path(logfile).absolute}'.")


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
        extra_flags:List[str]=None,
        output_dir:str='./',
        delete_scores:bool=True,
        delete_crash:bool=True,
        ) -> pd.DataFrame:
        """Runs Rosetta's relax protocol on the supplied .pdb file, returning the path to the relaxaed structure as 
        well as a dictionary with all characteristics found in the score.sc file. Function provides direct access to a 
        number of commandline options and the ability to add arbitrary commandline options at the end. 

        NOTE: The majority of crashes from this function are due to bad atom names.
         

        Args:
            infile: A str() with the path to the .pdb file to relax. 
            nstruct: Number of structures to create. 
            ignore_zero_occupancy: If relax should ignore atoms with zero occupancy. True by default.
            full_atom: If relax should do full atom relaxation. True by default.
            detect_disulf: If Rosetta should detect disulfide bonds. True by default.
            linmem_ig: Number of recent rotamers to store. 10 by default.
            constrain_relax_to_start_coords: If the backbone atoms should be constrained. True by default.
            coord_constrain_sidechains: If the sidechain heavy atoms should be constrained. True by default.
            ramp_constraints: If the constraints should be ramped during initial relaxation stage. True by default.
            prefix: str() with prefix for output file names. None and not used by default.
            overwrite: If results should be overwritten. True by default.
            extra_flags: A List[str] of extra flags to be added to the commandline. Empty by default. NOT CHECKED FOR CORRECTNESS. 
            output_dir: The output directory where the files will be saved. './' by default.
            delete_scores: Whether the score.sc file should be deleted after running. True by default.
            delete_crash: Whether the ROSETTA_CRASH.log file should be deleted after running. True by default.


        Returns:
            pandas DataFrame containing the results and energies of the relaxed structures. Description column contains
            full paths to relaxed files. 
        """
        #TODO(CJ):should also be able to take a Structure as input.
        fs.check_file_exists( infile )

        if Path(infile).suffix != '.pdb':
            _LOGGER.error(f"Expected input file format is .pdb. {infile} is an invalid entry. Exiting...")
            exit( 1 )

        fs.safe_rm( f'{output_dir}/score.sc')
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
        flags.append(f"-out:path:all {output_dir}")


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

        fs.safe_mkdir( output_dir )

        self.env_manager_.run_command(
            self.config_.RELAX,
            flags
        )

        df:pd.DataFrame = self.parse_score_file(f'{output_dir}/score.sc')

        df['description'] = df.apply( lambda row: f"{output_dir}/{row.description}.pdb", axis=1 )


        if delete_scores:
            self._delete_score_file( output_dir )

        if delete_crash:
            self._delete_crash_log( output_dir )


        return df
            


    def score(self,
            infile:str,
            ignore_zero_occupancy:bool=True,
            overwrite:bool=True,
            extra_flags:List[str]=None,
            output_dir:str='./',
            delete_scores:bool=True,
            delete_crash:bool=True,

        ) -> float:
        """Provides the total score in Rosetta Energy Units (REU) for a given structure. Uses default flags but can have behavior modified
        via supplied extra_flags. Returns the total score in REU.

        Arguments:
            infile: A str() with the path to the .pdb file to relax. 
            ignore_zero_occupancy: If relax should ignore atoms with zero occupancy. True by default.
            overwrite: If results should be overwritten. True by default.
            extra_flags: A List[str] of extra flags to be added to the commandline. Empty by default. NOT CHECKED FOR CORRECTNESS. 
            output_dir: The output directory where the files will be saved. './' by default.
            delete_scores: Whether the score.sc file should be deleted after running. True by default.
            delete_crash: Whether the ROSETTA_CRASH.log file should be deleted after running. True by default.

        Returns:
            Score of structure in file in REU.

        """
        fs.check_file_exists( infile )
        
        flags:List[str]=[
            f"-in:file:s '{infile}'",
             "-ignore_unrecognized_res",
        ]
        
        flags.append(f"-ignore_zero_occupancy {'true' if ignore_zero_occupancy else 'false'}")
        flags.append(f"-out:path:all {output_dir}")

        if overwrite:
            flags.append("-overwrite")

        if extra_flags:
            flags.extend(extra_flags)

        fs.safe_rm( f"{output_dir}/score.sc" )

        fs.safe_mkdir( output_dir )

        self.env_manager_.run_command(
            self.config_.SCORE,
            flags
        )


        df:pd.DataFrame=self.parse_score_file(f"{output_dir}/score.sc")

        if len(df) != 1:
            _LOGGER.error("Found more than one entry in score.sc file. Exiting...")
            exit( 1 )

        if delete_scores:
            self._delete_score_file( output_dir )

        if delete_crash:
            self._delete_crash_log( output_dir )

        return df.iloc[0].total_score


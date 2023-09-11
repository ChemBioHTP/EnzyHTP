"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation,
relaxation (minimization), scoring, ligand parameterization, and the ability to use RosettaScripts.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
import re
import shutil
from pathlib import Path
from copy import deepcopy
from xml.dom import minidom
import xml.etree.cElementTree as ET
from collections import namedtuple
from typing import Any, Dict, List, Tuple, Set, Union

import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em

from enzy_htp._config.rosetta_config import RosettaConfig, default_rosetta_config

from .base_interface import BaseInterface


class RosettaCst:
    """Class that reprsents a constraint set in the Rosetta molecular modelling package. Stores information about the residues
    and atoms involved. 
    TODO(CJ):

    """
    
    ALLOWED_CSTS: Set[str] = {"distanceAB", "angle_A", "angle_B", "torsion_A", "torsion_B", "torsion_AB"}
    """ """        

    def __init__(self, 
                parent,
                rname_1:str, 
                rnum_1:int, 
                ratoms_1:List[str], 
                rchain_1:str,
                rname_2:str,
                rnum_2:int,
                ratoms_2:List[str],
                rchain_2:str,
                constraints:List[Union[str,int]]):
        """"""
        self.parent_ = parent
        self.rname_1 = rname_1
        self.rnum_1 = rnum_1 
        self.ratoms_1 = ratoms_1 
        self.rchain_1 = rchain_1 
        self.rname_2 = rname_2 
        self.rnum_2 = rnum_2 
        self.ratoms_2 = ratoms_2 
        self.rchain_2 = rchain_2 
        self.constraints = constraints 

        #TODO(CJ): do some kind of validation for the constraints


    def parent(self):
        return self.parent_

    def set(self, key:str, value: Any) -> None:
        """Sets an attribute specified by the 'key' to the given 'value'. No checks are performed
        on the supplied value.

        Args:
            key: The str() name of the attribute to change.
            value: What you want to set the attribute to.

        Returns:
            Nothing.

        Raises:
            KeyError() if the supplied attribute does not exist. 
        """
        if key not in self.__dict__:
            raise KeyError(f"{key} is not in the constraint")

        self.__dict__[key] = value


    def evaluate(self, file:str) -> List[float]:
        """ 

        Args:
            file: Name of the file to analyze.

        Returns:
            
        """
        session = self.parent().parent().pymol.new_session()
        self.parent().parent().pymol.general_cmd(session,[
            ('delete','all'),
            ('load', file)
        ])

        sele1, sele2, sele3, sele4 = None, None, None, None

        differences:List[float] = list() 

        for cst in self.constraints:
            #print(cst)
            #TODO(CJ): do a tuple expansion here
            cst_type:str = cst[0]
            target = cst[1]
            tolerance = cst[2]
            #TODO(CJ): need to check if the angle is weird for this

            if cst_type.startswith('distance'):
                sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                dist:float = self.parent().parent().pymol.general_cmd(session, [('distance', None, sele1, sele2)] )[0]
                differences.append( abs(dist - target ) / tolerance) 
            elif cst_type.startswith('angle'):
                #TODO(CJ): need to address periodicity 
                if cst_type == 'angle_A':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                elif cst_type == 'angle_B':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                angle:float = self.parent().parent().pymol.general_cmd(session,[('angle', None, sele1, sele2, sele3)])[0]
                differences.append( abs(angle-target)/tolerance )
            elif cst_type.startswith('dihedral'):
                if cst_type == 'torsionA':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[2]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele3 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                elif cst_type == 'torsionAB':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[1]}"
                    sele2 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                elif cst_type == 'torsionB':
                    sele1 = f"chain {self.rchain_1} and resi {self.rnum_1} and name {self.ratoms_1[0]}"
                    sele2 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[0]}"
                    sele3 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[1]}"
                    sele4 = f"chain {self.rchain_2} and resi {self.rnum_2} and name {self.ratoms_2[2]}"
                assert False
                args.append(
                    ('dihedral', sele1, sele2, sele3, sele4 )
                )

        return differences

    def contains(self, chain:str, res_num:int) -> bool:
        """Does the constraint contain the residue at <chain>.<res_num>?"""
        if chain == self.rchain_1 and res_num == self.rnum_1:
            return True

        if chain == self.rchain_2 and res_num == self.rnum_2:
            return True

        return False

    def create_pdb_line(self, idx:int) -> str:
        """TODO(CJ): documentation"""
        return f"REMARK 666 MATCH TEMPLATE {self.rchain_1} {self.rname_1}  {self.rnum_1:>3} MATCH MOTIF {self.rchain_2} {self.rname_2}  {self.rnum_2:>3}  {idx:>3}  1"

    def create_cst_lines(self) -> List[str]:
        """TODO(CJ): documentation"""
        cst_content:List[str] = list()
        cst_content.append("CST::BEGIN")
        cst_content.append(
            f"   TEMPLATE::  ATOM_MAP: 1 atom_name: {' '.join(self.ratoms_1)}"
        )
        cst_content.append(
            f"   TEMPLATE::  ATOM_MAP: 1 residue3: {self.rname_1}")
        cst_content.append("")
        cst_content.append(
            f"   TEMPLATE::  ATOM_MAP: 2 atom_name: {' '.join(self.ratoms_2)}"
        )
        cst_content.append(
            f"   TEMPLATE::  ATOM_MAP: 2 residue3: {self.rname_2}")
        cst_content.append("")
        
        for ridx, rule in enumerate(self.constraints):
            if rule[0] == 'distanceAB':
                end = 0
            else:
                end = rule[4]
            cst_content.append(
                f"   CONSTRAINT::  {rule[0]:>10}: {rule[1]:6.2f} {rule[2]:6.2f} {rule[3]:6.2f} {end}"
            )
        
        cst_content.append("CST::END")
        
        return cst_content


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

    def _delete_score_file(self, working_dir: str = './') -> None:
        """Helper method that deletes the score.sc file in the specified directory.

        Args:
            working_dir: The directory to look for the score.sc file in. Defaults to './'

        Returns:
            Nothing.
        """
        score_file = working_dir + "/score.sc"
        fs.safe_rm(score_file)

    def _delete_crash_log(self) -> None:
        """Helper method that deletes the ROSETTA_CRASH.log file in the current working directory. ROSETTA_CRASH.log
        should always be in the current working directory.

        Args:
            Nothing. 

        Returns:
            Nothing.
        """
        fs.safe_rm('./ROSETTA_CRASH.log')

    def run_rosetta_scripts(self, opts: List[str], logfile: str = None) -> None:
        """Method that runs the rosettascripts executabl along with the supplied options. Optionally outputs
        the stdout to a specified logfile. Note that no sanitation is performed prior to running.

        Args:
            opts: a list() of str() to be run by the RosettaScripts executable.
            logfile: The file to output the stdout log to. Optional.

        Returns:
            Nothing,
        """

        if logfile:
            opts.extend([">", str(logfile)])

        self.env_manager_.run_command(self.config_.ROSETTA_SCRIPTS, opts)

        if logfile:
            _LOGGER.info(f"Saved RosettaScripts log to '{Path(logfile).absolute}'.")

    def parse_score_file(self, fname: str) -> pd.DataFrame:
        """Method that parses a score file into a Pandas Dataframe. Only encodes lines that begin with SCORE.

        Args:
            fname: Path to the score file. Will error if does not exist.

        Returns:
            A pandas dataframe containing the data in the supplied score file.
        """

        if not Path(fname).exists():
            _LOGGER.error(f"The suppliied file '{fname}' does not exist. Exiting...")
            exit(1)

        lines: List[str] = fs.lines_from_file(fname)

        lines = list(filter(lambda ll: ll.startswith('SCORE:'), lines))

        data = list(map(lambda ll: ll.split()[1:], lines))

        df = pd.DataFrame(data=data[1:], columns=data[0])

        column_names = list(df.columns)
        for cn in column_names:
            if cn == 'description':
                continue
            df[cn] = df[cn].astype('float')

        return df

    def parameterize_ligand(self, molfile: str, res_name: str, outdir: str = None, conformers: str = None) -> Tuple[str, str]:
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
        fs.check_file_exists(molfile)

        ALLOWED_FORMATS: List[str] = ".sdf .mol2 .mol".split()

        if Path(molfile).suffix not in ALLOWED_FORMATS:
            _LOGGER.error(
                f"The supplied file '{molfile}' is of an unsupported file type. Supported types are {', '.join(ALLOWED_FORMATS)}. Exiting..."
            )
            exit(1)
        #TODO(CJ): need to fix this. codes like 152 are valid  => only numbers
        #if not res_name.isupper() or not len(res_name) == 3 or not res_name.isalpha():
        #    _LOGGER.error(f"The supplied residue name '{res_name}' is invalid. It must be alphanumeric, capitalized, and have three characters. Exiting...")
        #    exit( 1 )
        #TODO(CJ): need to add some flags here -> potentially the charge one
        flags: List[str] = [self.config_.PARAMS_SCRIPT, f"{molfile}", f"--name={res_name}", "--clobber", "--keep-names"]

        self.env_manager_.run_command(self.config_.PY_2_7, flags)

        params_file: str = f"./{res_name}.params"
        pdb_file: str = f"./{res_name}_0001.pdb"

        if outdir:
            fs.safe_mkdir(outdir)

            outdir = Path(outdir)
            params_start = params_file
            pdb_start = pdb_file

            params_file = str(outdir / Path(params_start).name)
            pdb_file = str(outdir / Path(pdb_start).name)
            #TODO(CJ): use the safe_mv function here
            shutil.move(params_start, params_file)
            shutil.move(pdb_start, pdb_file)

        fs.check_file_exists(params_file)
        fs.check_file_exists(pdb_file)

        if conformers:
            #TODO(CJ): figure out exactly how I deal with the movment of conformer files
            self.add_conformers(res_name, params_file, pdb_file, conformers )

        return (params_file, pdb_file)


    def _fix_conformers(self, res_code:str, pdb_template:str, conformers_file:str) -> str:
        """TODO(CJ)"""
        session = self.parent().pymol.new_session()
        self.parent().pymol.general_cmd(session, [("load", pdb_template)])
        template_df:pd.DataFrame = self.parent().pymol.collect(session, "memory", "name elem x y z ID chain resn resi".split())
        
        original:str = self.parent().pymol.general_cmd(session, [
            ("delete", "all"), ("load", conformers_file), ("get_object_list",'all')
        ])[-1]

        #original = session.cmd.get_object_list()
        assert len(original) == 1
        split:List[str] = self.parent().pymol.general_cmd(session, [
            ("split_states","all"), ("get_object_list","all")
        ])[-1]
        print(split)
#        original = original[0]
#        content: List[str] = list()
#        #TODO(CJ): figure out a way to get rid of the warnings
#        #redirect_stdout()
        for oidx, oo in enumerate(split):
            if oo == original:
                continue
            df:pd.DataFrame = self.parent().pymol.collect(session, "memory", "name elem x y z ID chain resn resi".split(), sele=oo)
            assert len(df) == len(template_df), f"{len(df)} {len(template_df)}"
            for (tidx, trow), (idx, row) in zip(template_df.iterrows(),
                                                df.iterrows()):
                assert trow.elem == row.elem
                self.parent().pymol.general_cmd(session,[
                    ("alter", f"{oo} and ID {row.ID}", "name=" + trow["name"] ),
                    ("alter", f"{oo} and ID {row.ID}", f"chain='{trow.chain}'"),
                    ("alter", f"{oo} and ID {row.ID}", f"resn='{res_code}'"),
                ])


            temp_fname = f"state_{oidx}.pdb"
            self.parent().pymol.general_cmd(session,[
                    ("save", temp_name, oo)
            ])

            for ll in fs.lines_from_file(temp_fname):
                if ll.startswith('HETATM') or ll.startswith('ATOM'):
                    content.append(ll)

            content.append('TER')
            fs.safe_rm(temp_fname)

        content.pop()
        content.append('END')

        outfile = Path(conformers).with_suffix('.pdb')
        fs.write_lines(outfile, content)
        print(outfile)
        exit( 0 )
        return str(outfile)


    def add_conformers(self, res_code:str, param_file: str, pdb_template:str, conformers_file: str) -> None:
        """TODO(CJ)"""

        fs.check_file_exists(param_file)
        fs.check_file_exists(conformers_file)

        suffix:str = Path(conformers_file).suffix
        if suffix != ".pdb":
            if suffix in ".mol2 .mol .sdf".split():
                #TODO(CJ): this does not work
                conformers_file = self._fix_conformers(res_code, pdb_template, conformers_file)
                pass
            else:
                _LOGGER.error(f"The supplied file '{conformers_file}' is not in the .pdb and cannot be converted from .mol2, .mol, or .sdf. Exiting...")
                exit(1)

        content: List[str] = fs.lines_from_file(param_file)
        content.append(f"PDB_ROTAMERS {Path(conformers_file).absolute()}")
        fs.safe_rm(param_file)
        fs.write_lines(param_file, content)

    def relax(
        self,
        infile: str,
        nstruct: int,
        ignore_zero_occupancy: bool = True,
        full_atom: bool = True,
        detect_disulf: bool = True,
        linmem_ig: int = 10,
        constrain_relax_to_start_coords: bool = True,
        coord_constrain_sidechains: bool = True,
        ramp_constraints: bool = True,
        prefix: str = None,
        overwrite: bool = True,
        extra_flags: List[str] = None,
        output_dir: str = './',
        delete_scores: bool = True,
        delete_crash: bool = True,
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
        fs.check_file_exists(infile)

        if Path(infile).suffix != '.pdb':
            _LOGGER.error(f"Expected input file format is .pdb. {infile} is an invalid entry. Exiting...")
            exit(1)

        fs.safe_rm(f'{output_dir}/score.sc')
        #/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/relax.default.linuxgccrelease
        #-out:prefix $prefix
        #-out:file:scorefile ${prefix}.sc &
        flags: List[str] = [
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

        fs.safe_mkdir(output_dir)

        self.env_manager_.run_command(self.config_.RELAX, flags)

        df: pd.DataFrame = self.parse_score_file(f'{output_dir}/score.sc')

        df['description'] = df.apply(lambda row: f"{output_dir}/{row.description}.pdb", axis=1)

        if delete_scores:
            self._delete_score_file(output_dir)

        if delete_crash:
            self._delete_crash_log()

        return df

    def score(
        self,
        infile: str,
        ignore_zero_occupancy: bool = True,
        overwrite: bool = True,
        extra_flags: List[str] = None,
        output_dir: str = './',
        delete_scores: bool = True,
        delete_crash: bool = True,
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
        fs.check_file_exists(infile)

        flags: List[str] = [
            f"-in:file:s '{infile}'",
            "-ignore_unrecognized_res",
        ]

        flags.append(f"-ignore_zero_occupancy {'true' if ignore_zero_occupancy else 'false'}")
        flags.append(f"-out:path:all {output_dir}")

        if overwrite:
            flags.append("-overwrite")

        if extra_flags:
            flags.extend(extra_flags)

        fs.safe_rm(f"{output_dir}/score.sc")

        fs.safe_mkdir(output_dir)

        self.env_manager_.run_command(self.config_.SCORE, flags)

        df: pd.DataFrame = self.parse_score_file(f"{output_dir}/score.sc")

        if len(df) != 1:
            _LOGGER.error("Found more than one entry in score.sc file. Exiting...")
            exit(1)

        if delete_scores:
            self._delete_score_file(output_dir)

        if delete_crash:
            self._delete_crash_log()

        return df.iloc[0].total_score

    def relax_loops(
        self,
        infile: str,
        nstruct: int,
        ignore_zero_occupancy: bool = True,
        full_atom: bool = True,
        detect_disulf: bool = True,
        linmem_ig: int = 10,
        constrain_relax_to_start_coords: bool = True,
        coord_constrain_sidechains: bool = True,
        ramp_constraints: bool = True,
        prefix: str = None,
        overwrite: bool = True,
        extra_flags: List[str] = None,
        output_dir: str = './',
        delete_scores: bool = True,
        delete_crash: bool = True,
    ) -> pd.DataFrame:
        """

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
        fs.check_file_exists(infile)

        if Path(infile).suffix != '.pdb':
            _LOGGER.error(f"Expected input file format is .pdb. {infile} is an invalid entry. Exiting...")
            exit(1)

        fs.safe_rm(f'{output_dir}/score.sc')
        #/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/relax.default.linuxgccrelease
        #-out:prefix $prefix
        #-out:file:scorefile ${prefix}.sc &
        flags: List[str] = [
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

        fs.safe_mkdir(output_dir)

        self.env_manager_.run_command(self.config_.RELAX, flags)

        df: pd.DataFrame = self.parse_score_file(f'{output_dir}/score.sc')

        df['description'] = df.apply(lambda row: f"{output_dir}/{row.description}.pdb", axis=1)

        if delete_scores:
            self._delete_score_file(output_dir)

        if delete_crash:
            self._delete_crash_log()

        return df

    def write_script(self, fname: str, args: List[Dict]) -> str:
        """Writes an XML script to be used with RosettaScripts. Each element of the XML script is represented
        as a dict() within a list() of args. Note that each element dict() is required to have two keys, "parent"
        and "tag". "parent" refers to which element it should be nested under. In the case that there are multiple
        levels of nesting, they are joined by "." characters. The "tag" is the tag name of the element but note
        that there is NO checking for whether or not the included element tags are valid.
        
        Args:
            fname: The .xml file to save the script contents to.
            args: A list() of dict()'s where each is an element in the final .xml file. 
            
        Returns:
            The name of the script file.
        """

        def _find_node(elem: ET.Element, name: str) -> ET.Element:
            """Helper function that recursively finds the specified parent XML node. Assumes that supplied name
            str() is correctly delimited with '.' characters. DOES NOT check for correctness of nam.e

            Args:
                elem: The ET.Element to search within.
                name: The str() name to search for.

            Return:
                The XML node with the target tag name.
            """
            tks: List[str] = name.split('.', 1)

            target: str = tks[0]

            result: ET.Element = None
            if elem.tag == target:
                result = elem 
            else:
                for ee in elem:
                    if ee.tag == target:
                        result = ee
                        break
                else:
                    _LOGGER.error(f"There is no element with tag name '{target}' at this level. Exiting...")
                    exit(1)

            if len(tks) > 1:
                return _find_node(result, tks[1])
            else:
                return result

        root = ET.Element("ROSETTASCRIPTS")
        ET.SubElement(root, "RESIDUE_SELECTORS")
        ET.SubElement(root, "SCOREFXNS")
        ET.SubElement(root, "LIGAND_AREAS")
        ET.SubElement(root, "INTERFACE_BUILDERS")
        ET.SubElement(root, "MOVEMAP_BUILDERS")
        #ET.SubElement(root, "SCORINGGRIDS")
        ET.SubElement(root, "TASKOPERATIONS")
        ET.SubElement(root, "SIMPLE_METRICS")
        ET.SubElement(root, "FILTERS")
        ET.SubElement(root, "MOVERS")
        ET.SubElement(root, "PROTOCOLS")
        ET.SubElement(root, "OUTPUT")

        for arg in args:
            parent_name = arg.pop("parent", None)
            tag_name = arg.pop("tag", None)
            target_node = None

            bad: bool = False

            if not parent_name:
                #TODO(CJ): check if the parent name has an equal sign specifying the target node
                # more than what you would specify it otherwise
                _LOGGER.error("No parent name supplied in XML element dict()!")
                bad = True

            if bad:
                _LOGGER.error("Problems with XML elements detected. Exiting...")
                exit(1)
            
            child_nodes = arg.pop('child_nodes', list())
            #if arg.get('append_elements_only', False):
            #    _ = arg.pop('append_elements_only')
            #    target_node = _find_node(root, tag_name)
            #    for attrib, value in arg.items():
            #        target_node.set( attrib, value )

            #else:
            parent: ET.Element = _find_node(root, parent_name)

            if arg.get('append_elements_only', False):
                _ = arg.pop('append_elements_only')
                if tag_name == "SCORINGGRIDS":
                    parent: ET.Element = _find_node(root, 'ROSETTASCRIPTS')
                    target_node = ET.Element( tag_name )
                    parent.insert(0, target_node)
                    #target_node = ET.SubElement(parent[0], tag_name )
                else:
                    target_node = ET.SubElement(parent, tag_name )

                for attrib, value in arg.items():
                    target_node.set( attrib, value )
            else:
                target_node = ET.SubElement(parent, tag_name, attrib=arg)


            
            if child_nodes:
                for cn in child_nodes:
                    tag_name = cn.pop('tag', None)
                    _ = cn.pop('parent', None)
                    #TODO(CJ): make this recursive so it actually works for super nested things
                    _ = ET.SubElement(target_node, tag_name, attrib=cn)
            

        for rr in root:
            rr.text = "\n\t"

        xmlstr: str = minidom.parseString(ET.tostring(root)).toprettyxml()
        xml_content: List[str] = xmlstr.replace('<?xml version="1.0" ?>\n', '').splitlines()
        
        xml_content = list(filter(
            lambda ll: len(ll.strip()) > 0,
            xml_content
        ))

        fs.write_lines(fname, xml_content)

        return fname

    def loop_relax(
        self,
        infile: str,
        nstruct: int,
        ignore_zero_occupancy: bool = True,
        detect_disulf: bool = True,
        linmem_ig: int = 10,
        overwrite: bool = True,
        extra_flags: List[str] = None,
        output_dir: str = './',
        delete_scores: bool = True,
        delete_crash: bool = True,
    ) -> pd.DataFrame:
        """TODO
        Args:
            infile: A str() with the path to the .pdb file to relax. 
            nstruct: Number of structures to create. 
            ignore_zero_occupancy: If relax should ignore atoms with zero occupancy. True by default.
            detect_disulf: If Rosetta should detect disulfide bonds. True by default.
            linmem_ig: Number of recent rotamers to store. 10 by default.
            overwrite: If results should be overwritten. True by default.
            extra_flags: A List[str] of extra flags to be added to the commandline. Empty by default. NOT CHECKED FOR CORRECTNESS. 
            output_dir: The output directory where the files will be saved. './' by default.
            delete_scores: Whether the score.sc file should be deleted after running. True by default.
            delete_crash: Whether the ROSETTA_CRASH.log file should be deleted after running. True by default.

        Returns:
            pandas DataFrame containing the results and energies of the relaxed structures. Description column contains
            full paths to relaxed files. 
        """
        fs.check_file_exists(infile)

        if Path(infile).suffix != '.pdb':
            _LOGGER.error(f"Expected input file format is .pdb. {infile} is an invalid entry. Exiting...")
            exit(1)

        fs.safe_rm(f'{output_dir}/score.sc')
        flags: List[str] = [
            f"-in:file:s '{infile}'",
            f"-nstruct {nstruct}",
            f"-linmem_ig {linmem_ig}",
        ]

        df: pd.DataFrame = self.parent().pymol.collect('production.pdb', 'resi ss resn'.split(), sele='name CA')
        df['resi'] = df.resi.astype(int)

        ss = []
        for i, row in df.iterrows():
            if row.resn.upper() in "MG ZN HG".split():
                ss.append('M')
            else:
                ss.append(row.ss)

        df['ss'] = ss

        elements: List[Dict] = [
            {
                'parent': 'SCOREFXNS',
                'tag': 'ScoreFunction',
                'name': 'score_fxn',
                'weights': 'ref2015'
            },
            {
                'parent': 'MOVERS',
                'tag': 'FastRelax',
                'name': 'fast_relax',
                'scorefxn': 'score_fxn'
            },
            {
                'parent': 'MOVERS.FastRelax',
                'tag': 'MoveMap',
                'name': 'move_map'
            },
        ]

        temp: Dict[str, str] = deepcopy({
            'state': df.iloc[0].ss,
            'start': df.iloc[0].resi,
            'end': df.iloc[0].resi,
        })

        for i, row in df.iterrows():
            #TODO(CJ): need to ignore the non-amino acid stuff here
            if row.ss == temp['state']:
                temp['end'] = row.resi
            else:
                elements.append(
                    deepcopy({
                        'parent': 'MOVERS.FastRelax.MoveMap',
                        'tag': 'Span',
                        'begin': str(temp['start']),
                        'end': str(temp['end']),
                        'chi': 'true',
                        'bb': 'true' if temp['state'] == 'L' else 'false'
                    }))
                temp = deepcopy({
                    'state': row.ss,
                    'start': row.resi,
                    'end': row.resi,
                })

        elements.append(
            deepcopy({
                'parent': 'MOVERS.FastRelax.MoveMap',
                'tag': 'Span',
                'begin': str(temp['start']),
                'end': str(temp['end']),
                'chi': 'true',
                'bb': 'true' if temp['state'] == 'L' else 'false'
            }))

        elements.append({'parent': 'PROTOCOLS', 'tag': 'Add', 'mover_name': 'fast_relax'})

        fpath = Path(infile)
        xml_input: str = fpath.parent / "__temp.xml"
        xml_script = self.write_script(xml_input, elements)

        flags.extend(['-parser:protocol', str(xml_input.absolute())])

        flags.append(f"-ignore_zero_occupancy {'true' if ignore_zero_occupancy else 'false'}")
        flags.append(f"-out:path:all {output_dir}")

        if detect_disulf:
            flags.append("-in:detect_disulf")

        if overwrite:
            flags.append("-overwrite")

        if extra_flags:
            flags.extend(extra_flags)

        fs.safe_mkdir(output_dir)

        self.run_rosetta_scripts(flags)

        df: pd.DataFrame = self.parse_score_file(f'{output_dir}/score.sc')

        df['description'] = df.apply(lambda row: f"{output_dir}/{row.description}.pdb", axis=1)

        if delete_scores:
            self._delete_score_file(output_dir)

        if delete_crash:
            self._delete_crash_log()

        return df


    def _parse_cst(self, raw:str) -> RosettaCst:
        """ """
        var:Dict[str,Any] = {'constraints': []}
    
        tokens:List[str] = list(filter(len,re.split('[()]',raw)))
    
        if len(tokens) < 3:
            _LOGGER.info(
                f"There must be at least 3 blocks in an individual constraint. There are only {len(tokens)} in '{raw}'. Exiting..."
            )
            exit(1)
    
        for tidx, tk in enumerate(tokens):
    
            spl: List[str] = tk.split(',')
            if tidx < 2:
                var[f"rchain_{tidx+1}"] = spl[0]
                var[f"rnum_{tidx+1}"] = int(spl[1])
                var[f"rname_{tidx+1}"] = spl[2]
                var[f"ratoms_{tidx+1}"] = spl[3:]
            else:
                cst_type:str=spl[0]
                if cst_type not in RosettaCst.ALLOWED_CSTS:
                    _LOGGER.error(
                        f"The supplied constraint type {cst_type} is not supported. Allowed are: {', '.join(sorted(list(RosettaCst.ALLOWED_CSTS)))}. Exiting..."
                    )
                    exit(1)
                
                temp = [cst_type]
                
                for tt in spl[1:]:
                    if tt.find('.') == -1:
                        temp.append(int(tt))
                    else:
                        temp.append(float(tt))
                
                t_len:int=len(temp)
    
                if t_len < 5:
    
                    #TODO(CJ): add in some discussion of using a database here
                    # and mention that we are using default parameters
                    if cst_type == 'distanceAB':
                        temp.extend(
                            [2.00,0.25,100.00,0][t_len-1:]
                        )
                    elif cst_type in 'angle_A angle_B'.split():
                        temp.extend(
                            [180.0,5.0,100.0,360.0,1][t_len-1:]
                        )
    
                    else:
                        raise TypeError()
    
                var['constraints'].append(temp)
        
        return RosettaCst(
                            parent=self,
                            rname_1=var['rname_1'],
                            rnum_1=var['rnum_1'],
                            ratoms_1=var['ratoms_1'],
                            rchain_1=var['rchain_1'],
                            rname_2=var['rname_2'],
                            rnum_2=var['rnum_2'],
                            ratoms_2=var['ratoms_2'],
                            rchain_2=var['rchain_2'],
                            constraints=var['constraints']
                        )
    
    
    def csts_from_file(self, fname:str) -> List[RosettaCst]:
        """ """
        fs.check_file_exists( fname )
    
        raw:str = fs.content_from_file( fname ) 
    
        return self.csts_from_str( raw )
    
    
    def csts_from_str(self, raw: str) -> List[RosettaCst]:
        """ """
        raw = ''.join(raw.split())
        
        if raw.count('(') != raw.count(')'):
            _LOGGER.error(f"Unbalanced parantheses in raw constraint '{raw}'. Exiting...")
            exit( 1 )
    
        result = list()
    
        for token in raw.split('),('):
            if token[0] != '(':
                token = '(' + token
    
            if token[-1] != ')':
                token = ')' + token
    
            result.append(
                self._parse_cst( token )
            )
    
    
        return result        






"""Defines a RosettaInterface class that serves as a bridge for enzy_htp to utilize the Rosetta modelling
suite. Uses the RosettaConfig class found in enzy_htp/_config/rosetta_config.py. Supported operations include mutation,
relaxation (minimization), scoring, ligand parameterization, and the ability to use RosettaScripts.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-03-28
"""
import shutil
from pathlib import Path
from copy import deepcopy
from xml.dom import minidom
import xml.etree.cElementTree as ET
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

        flags: List[str] = [self.config_.PARAMS_SCRIPT, f"{molfile}", f"--name={res_name}", "--clobber"]

        self.env_manager_.run_command("python2.7", flags)

        params_file: str = f"./{res_name}.params"
        pdb_file: str = f"./{res_name}_0001.pdb"

        if outdir:
            fs.safe_mkdir(outdir)

            outdir = Path(outdir)
            params_start = params_file
            pdb_start = pdb_file

            params_file = str(outdir / Path(params_start).name)
            pdb_file = str(outdir / Path(pdb_start).name)

            shutil.move(params_start, params_file)
            shutil.move(pdb_start, pdb_file)

        fs.check_file_exists(params_file)
        fs.check_file_exists(pdb_file)

        if conformers:
            self.add_conformers(params_file, conformers)

        return (params_file, pdb_file)

    def add_conformers(self, param_file: str, conformers_file: str) -> None:
        """Adds conformers to the end of the .params file for a given ligand. Checks that files exist and
        exits if not.

        Args:
            param_file: The .params file as a str().
            conformers_file: The path to the conformers file. MUST BE IN .pdb FORMAT. 

        Returns:
            Nothing.
        """

        fs.check_file_exists(param_file)
        fs.check_file_exists(conformers_file)

        if not Path(conformers_file).suffix == ".pdb":
            _LOGGER.error(f"The supplied file '{conformers_file}' is not in the .pdb format. Exiting...")
            exit(1)

        content: List[str] = fs.lines_from_file(param_file)
        content.append(f"PDB_ROTAMERS {Path(conformers_file).absolute()}")
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
        ET.SubElement(root, "SCOREFXNS")
        ET.SubElement(root, "RESIDUE_SELECTORS")
        ET.SubElement(root, "TASKOPERATIONS")
        ET.SubElement(root, "SIMPLE_METRICS")
        ET.SubElement(root, "FILTERS")
        ET.SubElement(root, "MOVERS")
        ET.SubElement(root, "PROTOCOLS")
        ET.SubElement(root, "OUTPUT")

        for arg in args:
            parent_name = arg.pop("parent", None)
            tag_name = arg.pop("tag", None)

            bad: bool = False

            if not parent_name:
                _LOGGER.error("No parent name supplied in XML element dict()!")
                bad = True

            if not parent_name:
                _LOGGER.error("No tag name supplied in XML element dict()!")
                bad = True

            if bad:
                _LOGGER.error("Problems with XML elements detected. Exiting...")
                exit(1)

            parent: ET.Element = _find_node(root, parent_name)

            _ = ET.SubElement(parent, tag_name, attrib=arg)

        for rr in root:
            rr.text = "\n\t"

        xmlstr: str = minidom.parseString(ET.tostring(root)).toprettyxml()
        xml_content: List[str] = xmlstr.replace('<?xml version="1.0" ?>\n', '').splitlines()

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


    def mutate(self,
        molfile:str, 
        mutations:List[str,int,str]
        ) -> str:
        """Using Rosetta to mutate a protein

        Args:
            molfile:
            mutations:
    
        Returns:
            The mutated structure.
        """
        
        # validation 
        for row in mutations:
            if len(row) != 3:
                pass
        pass
        args:List[Dict] = list() 
        xml_file:str=f"{Path(mofile).parent}/__temp.xml"
        self.write_script(
            xml_file,
            args
        )
    

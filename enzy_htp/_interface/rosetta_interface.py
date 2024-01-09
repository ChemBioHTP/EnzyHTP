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
from plum import dispatch
from xml.dom import minidom
import xml.etree.cElementTree as ET
from collections import namedtuple
from typing import Any, Dict, List, Tuple, Set, Union

import numpy as np
import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em


from enzy_htp.structure import Structure, PDBParser, Mol2Parser, Ligand
from enzy_htp.structure.structure_constraint import StructureConstraint, ResiduePairConstraint

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

    def rename_atoms(self, stru: Structure) -> None:
        """Renames residues and atoms to be compatible with Rosetta naming and functions.
        
        Args:
            stru: The Structure() to perform renaming on.

        Returns:
            Nothing. 
        """
        nterm_mapper:Dict = {"H1":"1H", "H2":"2H", "H3":"3H"}
        his_mapper:Dict = {"HB2":"1HB", "HB3":"2HB"}
        _LOGGER.info("Beginning renaming...")
        changed_residues:int = 0
        changed_atoms:int = 0
        for res in stru.residues:
            if not res.is_canonical():
                continue
            
            if res.name in "HID HIS HIE".split():
                res.name = "HIS"
                changed_residues += 1
                for aa in res.atoms:
                    if aa.name in his_mapper:
                        aa.name = his_mapper[aa.name]
                        changed_atoms += 1

            for aa in res.atoms:
                if aa.name in nterm_mapper:
                    aa.name = nterm_mapper[aa.name]
                    changed_atoms += 1
        _LOGGER.info(f"Finished renaming! Changed {changed_residues} residues and {changed_atoms} atoms.")

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

        self.env_manager_.run_command(self.config_.ROSETTA_SCRIPTS, opts, quiet_fail=True)

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

    def parameterize_ligand(self, mol: Ligand, charge:int=None, work_dir:str=None) -> Tuple[str, str]:
        """Parameterizes the input ligand for use in the RosettaLigand protocol. Takes an input file with the ligand,
        as well as the name of the residue in PDB format (3 capitalized letters) as well as optionally where the output
        directory where the .params and .pdb files should be saved. The underlying script only supports .mol, .mol2, 
        and .sdf formats. Can also add conformers to end of .params file when conformers file is supplied. Function 
        exits on invalid inputs. 
        TODO(CJ): fix this documentation            
        Args:
            molfile: The name of the input file as a str().
            res_name: The all-capitalized, three letter string of the ligand in PDB format.
            outfile: Where the .params file should be saved. Optional.
            conformers: The conformers file for the given ligand. Optional.

        Returns:
            A tuple() with the layout of (.params file, .pdb file)            

        """
        if work_dir is None:
            work_dir = "./"

        res_name:str=mol.name
        molfile:str = f"{work_dir}/{res_name}.mol2"
        conformers:str=f"{work_dir}/{res_name}_conformers.pdb"
        params_file:str=f"./{res_name}.params"
        _parser = Mol2Parser()
        _parser.save_ligand(molfile, mol)
        flags: List[str] = [self.config_.PARAMS_SCRIPT, f"{molfile}", f"--name={res_name}", "--clobber", "--keep-names" ]
       
        indices = list()
        for aidx,aa in enumerate(mol.atoms):
            if aa.element == 'H':
                continue
            indices.append( aidx + 1 )
            
        if len(indices) <= 2:            
            flags.append( f"--nbr_atom={indices[0]}" )
        

        self.env_manager_.run_command(self.config_.PY_2_7, flags)
        fs.safe_rm(f"{res_name}_0001.pdb")
        
        params_file = fs.safe_mv(params_file, work_dir)
        params_content:List[str]=fs.lines_from_file(params_file)

        if charge is not None:
            params_content.append(
                f"NET_FORMAL_CHARGE {charge}"
            )

        _LOGGER.info(params_file)

        n_conformers:int=mol.n_conformers()
        if n_conformers > 1:
            conformer_file_content:List[str] = list()
            mol2_temp:str=f"{work_dir}/__temp_ligand_mol2.mol2"
            _LOGGER.info(f"Detected {n_conformers} in ligand {res_name}")
            _parser = Mol2Parser()
            session = self.parent().pymol.new_session()
            for cidx in range(1, n_conformers):
                conf = mol.get_ligand_conformer(cidx)
                fs.safe_rm(mol2_temp)
                _parser.save_ligand(mol2_temp, conf)
                pdb_conf:str=self.parent().pymol.convert(session, mol2_temp, new_ext='.pdb')
                conformer_file_content.extend(fs.lines_from_file(pdb_conf))
                fs.safe_rm(pdb_conf)

            conformer_file_content = list(filter(lambda ll: ll.startswith('HETATM') or ll.startswith('END'), conformer_file_content))

            fs.write_lines(conformers, conformer_file_content) 

            params_content.append(f"PDB_ROTAMERS {conformers}")

            fs.safe_rm(mol2_temp)

        fs.write_lines(params_file, params_content)
        
        return params_file 

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
                    target_node = ET.Element(tag_name)
                    parent.insert(0, target_node)
                    #target_node = ET.SubElement(parent[0], tag_name )
                else:
                    target_node = ET.SubElement(parent, tag_name)

                for attrib, value in arg.items():
                    target_node.set(attrib, value)
            else:
                target_node = ET.SubElement(parent, tag_name, attrib=arg)

            if child_nodes:
                for cn in child_nodes:
                    child_child_nodes = cn.pop('child_nodes', None)
                    tag_name = cn.pop('tag', None)
                    _ = cn.pop('parent', None)
                    #TODO(CJ): make this recursive so it actually works for super nested things
                    placed_child = ET.SubElement(target_node, tag_name, attrib=cn)

                    if child_child_nodes:
                        for ccn in child_child_nodes:
                            tag_name = ccn.pop('tag', None)
                            _ = ccn.pop('parent', None)
                            _ = ET.SubElement(placed_child, tag_name, attrib=ccn)

        for rr in root:
            rr.text = "\n\t"

        xmlstr: str = minidom.parseString(ET.tostring(root)).toprettyxml()
        xml_content: List[str] = xmlstr.replace('<?xml version="1.0" ?>\n', '').splitlines()

        xml_content = list(filter(lambda ll: len(ll.strip()) > 0, xml_content))

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


    def create_cst_pdb_line(self, cst:ResiduePairConstraint, idx: int) -> str:
        """Creates a str() PDB line in the appropriate format so that Rosetta can apply the constrained geometry
        described by the RosettaCst. This line goes in the corresponding PDB file.
        
        Args:
            idx: The index of the RosettaCst as an int().
            
        Returns:
            The PDB line corresponding to the RosettaCst.

        """
        return f"REMARK 666 MATCH TEMPLATE {cst.residue1.parent.name} {cst.residue1.name}  {cst.residue1.idx:>3} MATCH MOTIF {cst.residue2.parent.name} {cst.residue2.name}  {cst.residue2.idx:>3}  {idx:>3}  1"



    def create_cst_lines(self, cst:ResiduePairConstraint) -> List[str]:
        """Creates a List[str] which describes the constrained geometries in the required enzyme design format for Rosetta.
        These lines go into the corresponding .cst file."""
        cst_content: List[str] = list()
        cst_content.append("CST::BEGIN")
        cst_content.append(f"   TEMPLATE::  ATOM_MAP: 1 atom_name: {' '.join(map(lambda aa: aa.name, cst.residue1_atoms))}")
        cst_content.append(f"   TEMPLATE::  ATOM_MAP: 1 residue3: {cst.residue1.name}")
        cst_content.append("")
        cst_content.append(f"   TEMPLATE::  ATOM_MAP: 2 atom_name: {' '.join(map(lambda aa: aa.name, cst.residue2_atoms))}")
        cst_content.append(f"   TEMPLATE::  ATOM_MAP: 2 residue3: {cst.residue2.name}")
        cst_content.append("")

        for ridx, (rname, rule) in enumerate(cst.child_constraints):
            end = 0                
            if rule.is_angle_constraint() or rule.is_dihedral_constraint():
                end = 1
            cst_content.append(f"   CONSTRAINT::  {rname:>10}: {float(rule.target_value):6.2f} {float(rule.params['rosetta']['tolerance']):6.2f} {float(rule.params['rosetta']['penalty']):6.2f} {end}")

        cst_content.append("CST::END")

        return cst_content


    def write_constraint_file(self, stru:Structure, constraints:List[StructureConstraint], work_dir:str = None) -> str:
        #TODO(CJ): this!
        if work_dir is None:
            work_dir = "./"
    
        lines:List[str] = list()
        for cst in constraints:
            if cst.is_distance_constraint():
                assert False
                pass
            elif cst.is_angle_constraint():
                assert False
            elif cst.is_dihedral_constraint():
                assert False
            elif cst.is_residue_pair_constraint():
                for (cst_name, child_cst) in cst.child_constraints:
                    if child_cst.is_distance_constraint():
                        ridx_1:int=stru.absolute_index(child_cst.atoms[0].parent, indexed=1)
                        ridx_2:int=stru.absolute_index(child_cst.atoms[1].parent, indexed=1)
                        lines.append(
                            f"AtomPair {child_cst.atoms[0].name} {ridx_1} {child_cst.atoms[1].name} {ridx_2} LINEAR_PENALTY {child_cst.target_value:.2f} 0.00 {child_cst['rosetta']['tolerance']:.2f} {child_cst['rosetta']['penalty']:.2f}"
                        )
                    elif child_cst.is_angle_constraint():
                        ridx_1:int=stru.absolute_index(child_cst.atoms[0].parent, indexed=1)
                        ridx_2:int=stru.absolute_index(child_cst.atoms[1].parent, indexed=1)
                        ridx_3:int=stru.absolute_index(child_cst.atoms[2].parent, indexed=1)
                        lines.append(
                            f"Angle {child_cst.atoms[0].name} {ridx_1} {child_cst.atoms[1].name} {ridx_2} {child_cst.atoms[2].name} {ridx_3} LINEAR_PENALTY {np.radians(child_cst.target_value):.2f} 0.00 {np.radians(child_cst['rosetta']['tolerance']):.2f} {child_cst['rosetta']['penalty']/np.radians(1):.2f}"
                        )
                    else:
                        assert False

        fs.safe_mkdir(work_dir)
        fname:str = f"{work_dir}/constraints.cst"
        fs.write_lines(fname, lines )
        return fname 


    def integrate_enzdes_constraints(self, stru:Structure, constraints:List[StructureConstraint], work_dir:str=None) -> Tuple[str,str]:
        #TODO(CJ): update this

        if work_dir is None:
            work_dir = "./"

        fs.safe_mkdir(work_dir)

        _LOGGER.info("Beginning RosettaCst constraint integration...")
        parser = PDBParser()
        file_str = parser.get_file_str(stru, if_renumber=False, if_fix_atomname=False)
    
        pdb_content: List[str] = ["HEADER                                            xx-MMM-xx"]
        cst_content: List[str] = list()
        counter = 1
        for cidx, cst in enumerate(constraints):
            if cst.is_residue_pair_constraint():
                pdb_content.append(self.create_cst_pdb_line(cst, counter))
                cst_content.extend(self.create_cst_lines(cst))
                counter += 1
    
        pdb_file: str = f"{work_dir}/start.pdb"
        cst_file: str = f"{work_dir}/rdock.cst"
    
        if not Path(pdb_file).exists():
            fs.write_lines(pdb_file, pdb_content + file_str.splitlines())
    
        if not Path(cst_file).exists():
            fs.write_lines(cst_file, cst_content)
    
        _LOGGER.info("RosettaCst constraint integration successful! Relevant files:")
        _LOGGER.info(f"\t.pdb file: {Path(pdb_file).absolute()}")
        _LOGGER.info(f"\t.cst file: {Path(cst_file).absolute()}")
    
        return (pdb_file, cst_file)




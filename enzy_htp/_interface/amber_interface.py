"""Defines an AmberInterface class that serves as a bridge for enzy_htp to utilize AmberMD software. Uses the AmberConfig class
found in enzy_htp/molecular_mechanics/amber_config.py. Supported operations include mutation with tLEaP, MolDynParameterizer for MD,
and MolDynStep for modular MD steps that can be minimization, heating, constant pressure production, or constant pressure
equilibration.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Date: 2022-06-02
"""
import copy
import glob
from io import StringIO
import os
import re
import shutil
from pathlib import Path
from subprocess import CalledProcessError, CompletedProcess, SubprocessError
from typing import Generator, List, Tuple, Union, Dict, Any
from dataclasses import dataclass
import pandas as pd
from sympy import sympify
from collections import Iterable
from enzy_htp.structure.chain import Chain
from enzy_htp.structure.structure_region.api import create_region_from_residues, create_region_from_selection_pattern

from .base_interface import BaseInterface
from .handle_types import (
    MolDynParameterizer, MolDynParameter,
    MolDynStep,
    MolDynResult, MolDynResultEgg,)
from .ncaa_library import search_ncaa_parm_file
from .gaussian_interface import gaussian_interface

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import math_helper as mh
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.core.exception import AddPDBError, tLEaPError, AmberMDError
from enzy_htp.core.general import get_interval_str_from_list
from enzy_htp._config.amber_config import AmberConfig, default_amber_config
from enzy_htp.structure.structure_io import pdb_io, prmtop_io
from enzy_htp.structure.structure_constraint import (
    StructureConstraint,
    CartesianFreeze,
    merge_cartesian_freeze)
from enzy_htp.structure import StruSelection
from enzy_htp.structure import (
    Structure,
    Atom,
    Residue,
    Ligand,
    MetalUnit,
    ModifiedResidue,
    NonCanonicalBase,
    StructureEnsemble,
    PDBParser)
from enzy_htp import config as eh_config

class AmberParameter(MolDynParameter):
    """the Amber format MD parameter. Enforce the Amber parameter format
    in AmberMDStep.
    Attribute:
        inpcrd: the path of the .inpcrd input coordinate file
        prmtop: the path of the .prmtop parameter topology file"""

    def __init__(self, inpcrd_path: str, prmtop_path: str, ncaa_chrgspin_mapper: Union[Dict, None] = None):
        self._inpcrd = inpcrd_path
        self._prmtop = prmtop_path
        if ncaa_chrgspin_mapper is None:
            ncaa_chrgspin_mapper = {}
        self.ncaa_chrgspin_mapper = ncaa_chrgspin_mapper

    #region == property ==
    @property
    def engine(self) -> str:
        return "Amber"

    def get_solvated_structure(self) -> Structure:
        """get the solvated structure corresponding to the parameters"""
        pass # TODO do we really need this?

    @property
    def file_list(self) -> List[str]:
        """return a list of files that composes the parameter"""
        return [self._inpcrd, self._prmtop]

    @property
    def topology_file(self) -> str:
        """return the path of the topology file that composes the parameter"""
        return self._prmtop

    @property
    def topology_parser(self) -> callable:
        """return the parser object for topology_file"""
        return prmtop_io.PrmtopParser(self.ncaa_chrgspin_mapper).get_structure

    @property
    def input_coordinate_file(self) -> str:
        """return the path of the input_coordinate file that is the starting geom of the MD"""
        return self._inpcrd

    #endregion

    #region == checker ==
    def is_valid(self) -> bool:
        """check whether the parameter files represented by {self} is valid"""
        result = 1
        # file exist
        result *= Path(self._inpcrd).exists()
        result *= Path(self._prmtop).exists()
        # file size not zero
        result *= os.path.getsize(self._inpcrd) != 0 #TODO(CJ): use the function I already implemented for this
        result *= os.path.getsize(self._prmtop) != 0
        # TODO add upon need
        return bool(result)

    #endregion


class AmberParameterizer(MolDynParameterizer):
    """the MD parameterizer for Amber.
    Constructer:
        AmberInterface.build_md_parameterizer()
    Attributes: (configuration of the parameterization process)
        parent_interface
        engine
        charge_method
        resp_engine
        resp_lvl_of_theory
        ncaa_param_lib_path
        force_renew_ncaa_parameter
        ncaa_net_charge_engine
        ncaa_net_charge_ph
        solvate_box_type
        solvate_box_size
        gb_radii
        parameterizer_temp_dir
        additional_tleap_lines
        keep_tleap_in
    Methods:
        run()
        """

    RADII_MAPPER={
        1 : 'mbondi',
        2 : 'mbondi2',
        5 : 'mbondi2',
        7 : 'bondi',
        8 : 'mbondi3'}
    """{igb methods -> Radii keyword} mapping"""

    def __init__(
            self,
            interface,
            force_fields: List[str],
            charge_method: str,
            resp_engine: str,
            resp_lvl_of_theory: str,
            ncaa_param_lib_path: str,
            force_renew_ncaa_parameter: bool,
            ncaa_net_charge_engine: str,
            ncaa_net_charge_ph: float,
            solvate_box_type: str,
            solvate_box_size: float,
            gb_radii: int,
            parameterizer_temp_dir: str,
            additional_tleap_lines: List[str],
            keep_tleap_in: bool,) -> None:
        self._parent_interface: AmberInterface = interface
        self.force_fields = force_fields
        self.charge_method = charge_method
        self.resp_engine = resp_engine # TODO make this actually work
        self.resp_lvl_of_theory = resp_lvl_of_theory # TODO make this actually work
        self.ncaa_param_lib_path = ncaa_param_lib_path
        self.force_renew_ncaa_parameter = force_renew_ncaa_parameter # TODO make this actually work
        self.ncaa_net_charge_engine = ncaa_net_charge_engine # TODO make this actually work
        self.ncaa_net_charge_ph = ncaa_net_charge_ph # TODO make this actually work
        self.solvate_box_type = solvate_box_type
        self.solvate_box_size = solvate_box_size
        self.gb_radii = gb_radii
        self._parameterizer_temp_dir = parameterizer_temp_dir
        self.additional_tleap_lines = additional_tleap_lines
        self.keep_tleap_in = keep_tleap_in

    @property
    def engine(self) -> str:
        return "Amber"

    @property
    def parent_interface(self):
        return self._parent_interface

    @property
    def parameterizer_temp_dir(self) -> str:
        return self._parameterizer_temp_dir

    @parameterizer_temp_dir.setter
    def parameterizer_temp_dir(self, val: str) -> None:
        self._parameterizer_temp_dir = val

    def run(self, stru: Structure) -> AmberParameter:
        """the parameterizer convert stru to amber parameter (inpcrd, prmtop)"""
        # 0. set up paths
        result_inpcrd = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm.inpcrd")
        result_prmtop = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm.prmtop")
        temp_prmtop = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm_missing_pdb_info.prmtop")
        temp_ref_pdb = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm_ref_pdb.pdb")
        fs.safe_mkdir(self.parameterizer_temp_dir)

        # 1. check stru diversity
        diversity = stru.chemical_diversity
        _LOGGER.debug(f"diversity: {diversity}")

        # 2. extract and parameterize each special component
        ligand_parms = {}
        maa_parms = {}
        metalcenter_parms = {}
        gaff_type = self._check_gaff_type()
        if gaff_type == "GAFF2" and self.charge_method == "AM1BCC":
            _LOGGER.warning("found combinations of GAFF2 and AM1BCC is not recommend per Amber Manual! Consider GAFF2/ABCG2")
        if "ligand" in diversity:
            for lig in stru.ligands:
                if lig.name not in ligand_parms: # avoid repeating calculation
                    ligand_parms[lig.name] = self._parameterize_ligand(lig, gaff_type)
        if "modified_residue" in diversity:
            for maa in stru.modified_residue:
                if maa.name not in maa_parms:
                    maa_parms[maa.name] = self._parameterize_modified_res(maa, gaff_type)
        if "metalcenter" in diversity:
            _LOGGER.warning(
                "Support for paramization of metalcenter is not ready yet."
                " You can use AmberParameterizer().additional_tleap_lines"
                " & change relevant residue names to apply pre-made parameters"
                " for metal centers.")
            # all_used_residue_names = []
            # all_used_atom_types = []
            # for metal in stru.metalcenters:
            #     parms, used_residue_names, used_atom_types = self._parameterize_metalcenter(
            #         metal, ligand_parms, maa_parms)
            #     metalcenter_parms.update(parms)
            #     all_used_residue_names.extend(used_residue_names)
            #     all_used_atom_types.extend(used_atom_types)

        # 3. write the combining tleap.in
        tleap_content, temp_dry_pdb = self._write_combining_tleap_input(
                                            stru=stru,
                                            ligand_parms=ligand_parms,
                                            maa_parms=maa_parms,
                                            metalcenter_parms=metalcenter_parms,
                                            result_inpcrd=result_inpcrd,
                                            result_prmtop=temp_prmtop,
                                            additional_tleap_lines=self.additional_tleap_lines,)

        # 4. run tleap
        self.parent_interface.run_tleap(tleap_content, keep_in_file=self.keep_tleap_in)

        # 5. run add_pdb
        PDBParser().save_structure(temp_ref_pdb, stru, if_renumber=False, if_fix_atomname=False)
        self.parent_interface.run_add_pdb(temp_prmtop, result_prmtop, temp_ref_pdb)

        # 6. clean up
        fs.clean_temp_file_n_dir([
            temp_dry_pdb,
            temp_prmtop,
            temp_ref_pdb,
            self.parameterizer_temp_dir,
        ])

        return AmberParameter(result_inpcrd, result_prmtop, stru.ncaa_chrgspin_mapper)

    def _parameterize_ligand(self, lig: Ligand, gaff_type: str) -> Tuple[str, List[str]]:
        """parameterize ligand for AmberMD, use ncaa_param_lib_path for customized
        parameters. Multiplicity and charge information can be set in Ligand objects.
        TODO add a method for auto-assigning charge and spin
        Returns:
            mol_desc_path, [frcmod_path, ...]"""
        # san check
        if not gaff_type:
            _LOGGER.error("The structure contains non-canonical residue"
                          " (lig or maa) but GAFF/GAFF2 is not used!"
                          f" Check you force_fields. (current: {self.force_fields})")
            raise ValueError
        # init
        fs.safe_mkdir(self.ncaa_param_lib_path)
        target_method = f"{self.charge_method}-{gaff_type}"

        # 0. search parm lib
        mol_desc_path, frcmod_path_list = search_ncaa_parm_file(lig,
                                            target_method=target_method,
                                            ncaa_lib_path=self.ncaa_param_lib_path)

        if mol_desc_path:
            if frcmod_path_list:
                return mol_desc_path, frcmod_path_list
        else:
            # 1. generate mol2 if not found
            mol_desc_path = f"{self.ncaa_param_lib_path}/{lig.name}_{target_method}.mol2" # the search ensured no existing file named this
            self.parent_interface.antechamber_ncaa_to_moldesc(ncaa=lig,
                                                              out_path=mol_desc_path,
                                                              gaff_type=gaff_type,
                                                              charge_method=self.charge_method)

        # 2. run parmchk2 on the PDB
        # (this also runs when mol_desc exsit but not frcmod)
        frcmod_path = f"{self.ncaa_param_lib_path}/{lig.name}_{target_method}.frcmod"
        self.parent_interface.run_parmchk2(in_file=mol_desc_path,
                                           out_file=frcmod_path,
                                           gaff_type=gaff_type)

        return mol_desc_path, [frcmod_path]

    def _parameterize_modified_res(self, maa: ModifiedResidue, gaff_type: str) -> Tuple[str, List[str]]:
        """parameterize modified residues for AmberMD, use ncaa_param_lib_path for customized
        parameters. Multiplicity and charge information can be set in ModifiedResidue objects."""

        # san check
        if not gaff_type:
            _LOGGER.error("The structure contains non-canonical residue"
                          " (lig or maa) but GAFF/GAFF2 is not used!"
                          f" Check you force_fields. (current: {self.force_fields})")
            raise ValueError
    
        # init
        fs.safe_mkdir(self.ncaa_param_lib_path)
        target_method = f"{self.charge_method}-{gaff_type}"

        # 0. search parm lib - same as ligand
        mol_desc_path, frcmod_path_list = search_ncaa_parm_file(maa,
                                            target_method=target_method,
                                            ncaa_lib_path=self.ncaa_param_lib_path)

        if mol_desc_path:
            if frcmod_path_list:
                return mol_desc_path, frcmod_path_list
        else:
            # 1. generate ac if not found
            mol_desc_path = f"{self.ncaa_param_lib_path}/{maa.name}_{target_method}.ac" # the search ensured no existing file named this
            self.parent_interface.antechamber_ncaa_to_moldesc(ncaa=maa,
                                                              out_path=mol_desc_path,
                                                              gaff_type=gaff_type)
        # 2.1 fix the wrong atom type given by antechamber

        # 3. run prepgen on ac & mc get prepin
        prepin_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.prepin")

        # 4. run antechamber on prepin get mol2
        mol2_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.mol2")

        # 5. run parmchk2 twice on prepin get frcmod & frcmod2
        frcmod_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.frcmod")
        frcmod2_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.frcmod2")
        raise Exception("TODO")
        return mol2_path, [frcmod_path, frcmod2_path] # TODO make sure whether mol2 works or do we even need it?

    def _parameterize_metalcenter(self, metal: MetalUnit,
                                  ligand_parms: Dict[str, Tuple[str, List[str]]],
                                  maa_parms: Dict[str, Tuple[str, List[str]]],) -> List :
        """TODO maybe rewrote MCPB.py"""
        mcpb_path = f"{self.ncaa_param_lib_path}/mcpb"
        fs.safe_mkdir(mcpb_path)
        # 1. run MCPB.py 1st step - get involved residues from "small", new res name
        mcpbin_path = fs.get_valid_temp_name(
            f"{mcpb_path}/{metal.name}.mcpbin")
        # 2. search parm lib TODO: prob by coordinate

        # 3. run gaussian calculation

        # 4. run MCPB 2/3/4 step get mol2 files, frcmod file, new pdb, new atom types, bond lines
        new_residue_names = []
        new_atom_types = []

        # 5. store mol2, frcmod, bond lines in parm_dict
        parm_dict = {}
        raise Exception("TODO")
        return parm_dict, new_residue_names, new_atom_types

    def _write_combining_tleap_input(self,
            stru: Structure,
            ligand_parms: Dict,
            maa_parms: Dict,
            metalcenter_parms: Dict,
            result_inpcrd: str,
            result_prmtop: str,
            additional_tleap_lines: List[str],) -> Tuple[str, str]:
        """combine mol_desc and parm file of each noncanonical parts and make the content
        of the input file for tleap.
        Returns:
            (tleap.in content, temp_pdb_path)"""
        # init file path
        temp_pdb_path = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/tleap_combine.pdb")
        # make content
        lines = []
        for ff in self.force_fields:
            lines.append(f"source {ff}")

        # support for custom lines
        if additional_tleap_lines is not None:
            if (isinstance(additional_tleap_lines, Iterable)
                and not isinstance(additional_tleap_lines, str)
                and isinstance(additional_tleap_lines[0], str)
                ):
                lines.extend(additional_tleap_lines)
            else:
                _LOGGER.error(f"`additional_tleap_line` needs to be a list of str. current: {repr(additional_tleap_lines)}")
                raise TypeError

        # NCAA parts

        # # add atom types TODO
        # for xxx, new_atom_types in metalcenter_parms.values():
        #     pass

        # ligand
        for ncaa_name, (mol_desc_path, frcmod_path_list) in ligand_parms.items():
            lines.extend(
                self._make_ncaa_tleap_lines(ncaa_name, mol_desc_path, frcmod_path_list))
        # MAA
        for ncaa_name, (mol_desc_path, frcmod_path_list) in maa_parms.items():
            lines.extend(
                self._make_ncaa_tleap_lines(ncaa_name, mol_desc_path, frcmod_path_list))
            
        # PB/GBRadii setting
        if self.gb_radii:
            lines.append(f"set default PBRadii {self.RADII_MAPPER[self.gb_radii]}")

        # load PDB
        pdb_io.PDBParser().save_structure(outfile=temp_pdb_path,
                                          stru=stru,
                                          if_renumber=False)
        lines.extend([
            f"a = loadpdb {temp_pdb_path}",
            "center a",
            "addions a Na+ 0",
            "addions a Cl- 0",
            f"solvate{self.solvate_box_type} a TIP3PBOX {self.solvate_box_size}",
            f"saveamberparm a {result_prmtop} {result_inpcrd}",
            "quit",
            "",
        ])

        result = "\n".join(lines)

        return result, temp_pdb_path

    def _make_ncaa_tleap_lines(self,
            ncaa_name: str,
            mol_desc_path: str,
            frcmod_path_list: List[str]) -> List[str]:
        """make ncaa moldesc and params lines in tleap.in"""
        result = []
        # params
        for frcmod in frcmod_path_list:
            result.append(f"loadAmberParams {frcmod}")
        # mol desc
        if fs.get_file_ext(mol_desc_path) in [".prepin", ".prepi"]:
            result.append(f"loadAmberPrep {mol_desc_path}")
        elif fs.get_file_ext(mol_desc_path) in [".mol2"]:
            result.append(f"{ncaa_name} = loadmol2 {mol_desc_path}")
        else:
            _LOGGER.error(f"Got file with wrong extension {mol_desc_path}")
            raise ValueError
        return result

    def _check_gaff_type(self) -> str:
        """check the GAFF type used for parameterization.
        Return None if non of them are used."""
        for ff in self.force_fields:
            pattern_match = re.search(r"[Gg][Aa][Ff][Ff]2?", ff)
            if pattern_match:
                return pattern_match.group().upper()
        return None

@dataclass
class AmberMDResultEgg(MolDynResultEgg):
    """This class defines the data format of the MolDynResultEgg in Amber's case"""
    traj_path: str
    traj_log_path: str
    rst_path: str
    prmtop_path: str
    parent_job: ClusterJob


class AmberMDStep(MolDynStep):
    """the modular MD step of Amber.
    Constructer:
        AmberInterface.build_md_step()
    Attributes: (necessary information of the modular MD step)
        engine
        parent_interface
        length
        timestep
        minimize
        temperature
        thermostat
        constrain
        restart
        core_type
        cluster_job_config
        if_report
        record_period
        work_dir
    Methods:
        if_report()
        make_job()
        translate()
        try_merge_job()"""
    def __init__(self, 
                 interface,
                 name: str, 
                 length: float,
                 timestep: float,
                 minimize: bool,
                 temperature: Union[float, List[Tuple[float]]],
                 thermostat: str,
                 pressure_scaling: str,
                 constrain: List[StructureConstraint],
                 restart: bool,
                 core_type: str,
                 cluster_job_config: Dict,
                 if_report: bool,
                 record_period: float,
                 keep_in_file: bool,
                 work_dir: str,) -> None:
        self._parent_interface: AmberInterface = interface
        self.name = name
        self._length = length
        self._timestep = timestep
        self.minimize = minimize
        self._temperature = temperature
        self._thermostat = thermostat
        self.pressure_scaling = pressure_scaling
        self._constrain = constrain
        self.restart = restart
        self.core_type = core_type
        self.cluster_job_config = cluster_job_config
        self._if_report = if_report
        self.record_period = record_period
        self.keep_in_file = keep_in_file
        self._work_dir = work_dir

    # region == property == 
    @property
    def engine(self) -> str:
        return "Amber"

    @property
    def parent_interface(self):
        return self._parent_interface

    @property
    def if_report(self) -> bool:
        """whether the step reports the output"""
        return self._if_report

    @if_report.setter
    def if_report(self, val: bool):
        """setter for if_report"""
        self._if_report = val

    @property
    def work_dir(self) -> bool:
        """the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, val: bool):
        """setter for work_dir"""
        self._work_dir = val

    @property
    def constrain(self) -> List[StructureConstraint]:
        return self._constrain

    @constrain.setter
    def constrain(self, val: List[StructureConstraint]):
        self._constrain = val

    @property
    def length(self) -> float:
        return self._length

    @length.setter
    def length(self, val: float):
        self._length = val

    @property
    def temperature(self) -> float:
        return self._temperature

    @temperature.setter
    def temperature(self, val: float):
        self._temperature = val

    @property
    def thermostat(self) -> str:
        return self._thermostat

    @thermostat.setter
    def thermostat(self, val: str):
        self._thermostat = val

    @property
    def timestep(self) -> float:
        return self._timestep

    @timestep.setter
    def timestep(self, val: float):
        self._timestep = val

    @property
    def md_config_dict(self) -> Dict:
        """get md configuration related attributes in a dictionary form"""
        return {
            "name" : self.name,
            "length" : self.length,
            "timestep" : self.timestep,
            "minimize" : self.minimize,
            "temperature" : self.temperature,
            "thermostat" : self.thermostat,
            "pressure_scaling" : self.pressure_scaling,
            "constrain" : self.constrain,
            "restart" : self.restart,
            "if_report" : self.if_report,
            "record_period" : self.record_period,
            "mdstep_dir" : self.work_dir,
        }
    # endregion

    # methods
    def make_job(self, input_data: Union[AmberParameter, AmberMDResultEgg],
                path_rel_to: str = None,) -> Tuple[ClusterJob, MolDynResultEgg]:
        """the method that make a ClusterJob that runs the step
        Args:
            input_data: the input data of the step
            path_rel_to: 
                if this argument is specified, the file path inside of the .in file
                will be written as the relative path to {path_rel_to}. However, the
                path stored in python variables are still absolute or relative to cwd
                so that they can be valid."""
        # 1. parse input
        if isinstance(input_data, AmberParameter): # build from parameter
            coord = input_data.input_coordinate_file
            prmtop = input_data.topology_file
        elif isinstance(input_data, AmberMDResultEgg): # build form previous step
            coord = input_data.rst_path
            prmtop = input_data.prmtop_path
        else:
            _LOGGER.error("only allow AmberParameter or AmberMDResultEgg as `input_data`")
            raise TypeError

        # 2. make .in file
        fs.safe_mkdir(self.work_dir)
        temp_mdin_file = self._make_mdin_file(path_rel_to=path_rel_to)

        # 3. make cmd
        md_step_cmd, traj_path, mdout_path, mdrst_path = self._make_md_cmd(
            temp_mdin_file,
            prmtop,
            coord,
            path_rel_to=path_rel_to)

        # 4. assemble ClusterJob
        cluster = self.cluster_job_config["cluster"]
        res_keywords = self.cluster_job_config["res_keywords"]
        env_settings = cluster.AMBER_ENV[self.core_type.upper()]
        sub_script_path = fs.get_valid_temp_name(f"{self.work_dir}/submit_{self.name}.cmd")
        job = ClusterJob.config_job(
            commands = md_step_cmd,
            cluster = cluster,
            env_settings = env_settings,
            res_keywords = res_keywords,
            sub_dir = "./", # because path are relative
            sub_script_path = sub_script_path,
        )
        job.mimo = { # temp solution before having a 2nd MD engine
            "commands" : [md_step_cmd,],
            "env_settings" : env_settings,
            "res_keywords" : res_keywords,
            "sub_dir" : "./",
            "md_names" : [self.name],
            "work_dir" : self.work_dir,
            "temp_mdin": [temp_mdin_file],
        }

        # 5. make result egg
        result_egg = AmberMDResultEgg(
            traj_path=traj_path,
            traj_log_path=mdout_path,
            rst_path=mdrst_path,
            prmtop_path=prmtop,
            parent_job=job,
        )

        return (job, result_egg)

    def run(self, input_data: Union[MolDynParameter, MolDynResult]) -> MolDynResult:
        """the method that runs the step"""
        # 1. parse input
        if isinstance(input_data, AmberParameter): # build from parameter
            coord = input_data.input_coordinate_file
            prmtop = input_data.topology_file
        elif isinstance(input_data, MolDynResult) and input_data.source == "amber": # build form previous step
            coord = input_data.last_frame_file
            prmtop = input_data.last_frame_parser.prmtop_file
        else:
            _LOGGER.error("only allow AmberParameter or MolDynResult(source='amber') as `input_data`")
            raise TypeError

        # 2. make .in file
        fs.safe_mkdir(self.work_dir)
        temp_mdin_file = self._make_mdin_file()

        # TODO do we need to clear existing rst/nc files? to address issue #95? unit test it

        # TODO running sander/pmemd should go inside AmberInterface. Apply when another function needs it.
        # 3. make cmd
        md_step_cmd, traj_path, mdout_path, mdrst_path = self._make_md_cmd(
            temp_mdin_file,
            prmtop,
            coord)
        if not self.if_report:
            traj_path = None

        # 4. run cmd
        md_cmd_exe = md_step_cmd.split(" ")[0]
        md_cmd_args = md_step_cmd.split(" ")[1:]
        parent_interface: AmberInterface = self.parent_interface
        this_md_run = parent_interface.env_manager_.run_command(
            exe=md_cmd_exe,
            args=md_cmd_args,
        )
        self.check_md_error(traj_path, mdout_path, prmtop, this_md_run)

        # 5. make result
        traj_parser = AmberNCParser(
            prmtop_file=prmtop,
            interface=self.parent_interface).get_coordinates
        traj_log_parser = self.parent_interface.read_from_mdout
        last_frame_parser = AmberRSTParser(
            prmtop_file=prmtop,
            interface=self.parent_interface)

        result = MolDynResult(
            traj_file = traj_path,
            traj_parser = traj_parser, 
            traj_log_file = mdout_path,
            traj_log_parser = traj_log_parser, 
            last_frame_file = mdrst_path,
            last_frame_parser = last_frame_parser,
            source="amber",
        )

        # 6. clean up
        clean_up_target = ["./mdinfo"]
        if not self.keep_in_file:
            clean_up_target.append(temp_mdin_file)
        fs.clean_temp_file_n_dir(clean_up_target)

        return result


    def _make_mdin_file(self, path_rel_to: str = None,) -> str:
        """make a temporary mdin file.
        *path* is based on self.work_dir and self.name.
        If the file exists, will change the filename to {self.name}_{index}.in
        The index start from 1.
        *content* is based on attributes of the instance.
        
        Args:
            path_rel_to: 
                if this argument is specified, the file path inside of the .in file
                will be written as the relative path to {path_rel_to}. However, the
                path stored in python variables are still absolute or relative to cwd
                so that they can be valid.

        Return: the path of the temp mdin file."""
        # path
        temp_mdin_file_path = fs.get_valid_temp_name(f"{self.work_dir}/{self.name}.in")
        # content
        self.parent_interface.write_to_mdin(self.md_config_dict, temp_mdin_file_path, path_rel_to=path_rel_to)
        return temp_mdin_file_path

    def _make_md_cmd(self, temp_mdin_file, prmtop, coord, path_rel_to: str = None,) -> str:
        """compile the sander/pmemd cmd from config
        path_rel_to: 
            if this argument is specified, the file path inside of the .in file
            will be written as the relative path to {path_rel_to}. However, the
            path stored in python variables are still absolute or relative to cwd
            so that they can be valid."""
        num_cores = self.cluster_job_config["res_keywords"]["node_cores"]
        executable = self.parent_interface.get_md_executable(self.core_type, num_cores)
        mdout_path = f"{self.work_dir}/{self.name}.out"
        mdrst_path = f"{self.work_dir}/{self.name}.rst"
        traj_path = f"{self.work_dir}/{self.name}.nc"
        if path_rel_to:
            md_step_cmd= (
                f"{executable} "
                f"-O " # overwrite; note that AmberMDStep() with same name will overwrite each other
                f"-i {fs.relative_path_of(temp_mdin_file, path_rel_to)} " # mdin file
                f"-o {fs.relative_path_of(mdout_path, path_rel_to)} " # mdout file path
                f"-p {fs.relative_path_of(prmtop, path_rel_to)} " # prmtop path
                f"-c {fs.relative_path_of(coord, path_rel_to)} " # init coord
                f"-r {fs.relative_path_of(mdrst_path, path_rel_to)} " # last frame coord
                f"-ref {fs.relative_path_of(coord, path_rel_to)} " # reference coord used only when nmropt=1
                f"-x {fs.relative_path_of(traj_path, path_rel_to)} " # the output md trajectory
            )
        else:
            md_step_cmd= (
                f"{executable} "
                f"-O " # overwrite; note that AmberMDStep() with same name will overwrite each other
                f"-i {temp_mdin_file} " # mdin file
                f"-o {mdout_path} " # mdout file path
                f"-p {prmtop} " # prmtop path
                f"-c {coord} " # init coord
                f"-r {mdrst_path} " # last frame coord
                f"-ref {coord} " # reference coord used only when nmropt=1
                f"-x {traj_path} " # the output md trajectory
            )
        return md_step_cmd, traj_path, mdout_path, mdrst_path

    def translate(self, result_egg: AmberMDResultEgg) -> MolDynResult:
        """the method convert engine specific results to general output.
        will also clean up temp files."""
        traj_file = result_egg.traj_path
        traj_parser = AmberNCParser(
            prmtop_file=result_egg.prmtop_path,
            interface=self.parent_interface).get_coordinates
        traj_log_file = result_egg.traj_log_path
        traj_log_parser = self.parent_interface.read_from_mdout
        last_frame_file = result_egg.rst_path
        last_frame_parser = AmberRSTParser(
            prmtop_file=result_egg.prmtop_path,
            interface=self.parent_interface)
        
        # error check
        self.check_md_error(traj_file, traj_log_file, result_egg.prmtop_path, result_egg.parent_job)

        # clean up
        clean_up_target = ["./mdinfo"]
        if not self.keep_in_file:
            clean_up_target.extend(result_egg.parent_job.mimo["temp_mdin"])
            # TODO also clean up .rs file but how to obtain it?
            # parse mdin file?
        fs.clean_temp_file_n_dir(clean_up_target)

        return MolDynResult(
            traj_file = traj_file,
            traj_parser = traj_parser, 
            traj_log_file = traj_log_file,
            traj_log_parser = traj_log_parser, 
            last_frame_file = last_frame_file,
            last_frame_parser = last_frame_parser,
            source="amber",
        )

    def check_md_error( self,
                        traj: str,
                        traj_log: str,
                        prmtop: str,
                        stdstream_source: Union[ClusterJob,
                                                CompletedProcess,
                                                SubprocessError]):
        """a check for whether an error occurs in MD is needed everytime before generating
        a MolDynResult.
        Possible Amber error info places:
        1. stdout/stderr
        2. mdout file
        3. more analysis of comm problem from the traj (TODO) 
        The ultimate goal is to summarize each error type and give suggestions.
        Notes: designed here so that it can also checks based on config of the current step"""
        # only when MD is not complete, there is an error
        if not self.parent_interface.is_md_completed(traj_log):
            # collect error info
            error_info_list = []
            # 1. stdout stderr
            # error types: Fall of bus, periodic box has changed too much, illegel mem
            if isinstance(stdstream_source, ClusterJob):
                with open(stdstream_source.job_cluster_log) as f:
                    stderr_stdout = f.read()
                    error_info_list.append(f"stdout/stderr(from job log):{os.linesep*2}{stderr_stdout}")
            elif isinstance(stdstream_source, CompletedProcess):
                error_info_list.append(f"stdout:{os.linesep*2}{stdstream_source.stdout}")
                error_info_list.append(f"stderr:{os.linesep*2}{stdstream_source.stderr}")
            elif isinstance(stdstream_source, SubprocessError):
                error_info_list.append(f"stdout:{os.linesep*2}{stdstream_source.stdout}")
                error_info_list.append(f"stderr:{os.linesep*2}{stdstream_source.stderr}")
            else:
                _LOGGER.error("Only allow ClusterJob, CompletedProcess, SubprocessError as input types for `stdstream_source`")
                raise TypeError
            # 2. mdout file
            # have no experience that error pops here yet. Add when observed
            # 3. traj analysis 
            # TODO
            # Note that if_report = 0 will make traj_path = None
            # TODO give suggestion for each detected types of errors.

            _LOGGER.error(f"Amber MD didn't finish normally.{os.linesep}{os.linesep.join(error_info_list)}")
            raise AmberMDError(error_info_list)
    
    @classmethod
    def try_merge_jobs(cls, job_list: List[ClusterJob]) -> List[ClusterJob]:
        """the classmethod that merge a list of sequential jobs from MolDynStep to fewer jobs.
        Also check for jobs that overwrite each other.
        Jobs are mergable when they have identical environment and resource
        and are sequential in steps. The merge specifically merges the command of each job.
        Returns: a list of reduced jobs.
        TODO: This is not a good design.
            This method could be generalized in ClusterJob 
            but need to develop parser of sub_script_str
            we should probably do this when adding the 2nd MD engine."""
        result = [job_list[0]]
        for job in job_list[1:]:
            merged_job = result.pop() # expose the job for merge
            # check names
            if not set(job.mimo["md_names"]).isdisjoint(set(merged_job.mimo["md_names"])):
                _LOGGER.warning(f"Found md steps with same names! ({set(job.mimo['md_names'])} & {set(merged_job.mimo['md_names'])})"
                                " Their results will overwrite each other."
                                "We sugggest to change the name and make them unique if you care about their result.")
            # try merge
            mergable = (
                job.mimo["env_settings"],
                job.mimo["res_keywords"],
                job.mimo["sub_dir"],                
            )
            exposed_mergable = (
                merged_job.mimo["env_settings"],
                merged_job.mimo["res_keywords"],
                merged_job.mimo["sub_dir"],                
            )
            if mergable == exposed_mergable:
                merged_cmds = merged_job.mimo["commands"] + job.mimo["commands"]
                env_settings = merged_job.mimo["env_settings"]
                res_keywords = merged_job.mimo["res_keywords"]
                sub_dir = merged_job.sub_dir
                md_names = merged_job.mimo["md_names"] + job.mimo["md_names"]
                work_dir = merged_job.mimo["work_dir"]
                merged_mdin = merged_job.mimo["temp_mdin"] + job.mimo["temp_mdin"]
                sub_script_path = fs.get_valid_temp_name(f"{work_dir}/submit_{'_'.join(md_names)}.cmd")
                # update merged job
                merged_job = ClusterJob.config_job(
                    commands = merged_cmds,
                    cluster = job.cluster,
                    env_settings = env_settings,
                    res_keywords = res_keywords,
                    sub_dir = sub_dir,
                    sub_script_path = sub_script_path,
                )
                merged_job.mimo = {
                    "commands" : merged_cmds,
                    "env_settings" : env_settings,
                    "res_keywords" : res_keywords,
                    "sub_dir" : sub_dir,
                    "md_names" : md_names,
                    "work_dir" : work_dir,
                    "temp_mdin" : merged_mdin,
                }
                result.append(merged_job) # add the job back
            else:
                result.append(merged_job)
                result.append(job) # add unmergable

        return result


class AmberInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize AmberMD software.
    Main supported operations:
    MolDynStep: a modular MD step.
    MolDynParameterizer: parameterizer for MD.
    run_tleap: tleap interface and a series of tleap based operations like tleap_clean_up_stru.
    TODO: add more
    * Users should use this class as the only way to interact with any functionality
      in Amber or associated tools like tleap or cpptraj.

    Attributes:
        config_	: The AmberConfig() class which provides settings for both running Amber and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() that interface with the shell and ensure all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    AMBER_FILE_FORMAT_MAPPER = {
        ".prepin" : "prepi",
    }
    """the mapper for {ext:format} for file formats defined/used by amber
     that are different than the extension TODO add upon need"""

    SUPPORTED_CHARGE_METHOD_MAPPER = {
        "AM1BCC" : "bcc", "bcc" : "bcc",
        "RESP" : "resp", "resp" : "resp",
        "rc" : "rc",
    }
    """dictionary that maps keywords to charge method that current supported in run_antechamber()"""

    MD_THERMOSTAT_MAPPER = {
        "constant_energy" : 0,
        "berendsen" : 1,
        "anderson" : 2,
        "langevin" : 3,
        "oin" : 9,
        "sin-respa" : 10,
        "bussi" : 11,
    }
    """The mapper for mapping thermostat names to Amber keyword numbers"""

    MD_PRESSURE_SCALING_MAPPER = {
        "none" : 0,
        "isotropic" : 1,
        "anisotropic" : 2,
        "semiisotropic" : 3,
    }
    """The mapper for mapping pressure scaling names to Amber keyword numbers"""

    MD_RESTART_MAPPER = {True: (5, 1), False: (1, 0)}
    """The mapper for mapping restart - (ntx, irest) value"""

    MD_TIMESTEP_SHAKE_MAPPER = {0.000002: (2, 2), 0.000001: (1, 1)}
    """The mapper for mapping timestep - (ntc, ntf) value"""

    def __init__(self, parent, config: AmberConfig = None) -> None:
        """Simplistic constructor that optionally takes an AmberConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, default_amber_config)

    # == interface general ==
    def display_config(self) -> None:
        """Prints all settings for the object's AmberConfig() inteface to stdout using AmberConfig.display()."""
        self.config_.display()

    # == general amber app interface ==

    # -- file format --
    def get_file_format(self, fname: str) -> str:
        """determine the file type give file path"""
        ext = fs.get_file_ext(fname)
        if "frcmod" in ext:
            return "frcmod"
        return self.AMBER_FILE_FORMAT_MAPPER.get(ext, ext[1:])

    def convert_traj_to_nc(self, traj_path: str, out_path: str, topology_path: str = str()) -> None:
        """Convert the given trajectory file to the Amber `.nc` format file in the out_path.
        
        Args:
            traj_path (str): Path to the input trajectory file.
            out_path (str): Path to the output Amber `.nc` format trajectory file.
            topology_path (str, optional): Path to the topology file. Default empty.
        """
        in_format = self.get_file_format(traj_path)

        if in_format == "nc":
            fs.safe_cp(traj_path, out_path)
        elif in_format == "mdcrd":
            if (not topology_path or not fs.is_path_exist(topology_path)):
                _LOGGER.error(f"Topology filepath ({topology_path if len(topology_path) else 'N/A'}) doesn't exist.")
                raise FileNotFoundError
            contents = [
                f"parm {topology_path}",
                f"trajin {traj_path}",
                f"trajout {out_path}",
                "run",
                "quit"
            ]
            contents = "\n".join(contents)
            self.run_cpptraj(contents)
        else:
            # TODO for other formats. We probably need a central format to reduce the engineering overhead
            _LOGGER.error(f"found unsupported file format: {in_format} ({traj_path})")
            raise ValueError

    def convert_top_to_prmtop(self, fpath: str, out_path: str) -> None:
        """convert the given topology file to the Amber .prmtop file in the out_path"""
        in_format = self.get_file_format(fpath)

        if in_format == "prmtop":
            fs.safe_cp(fpath, out_path)
        else:
            _LOGGER.error(f"found unsupported file format: {in_format} ({fpath})")
            raise ValueError

    # -- tleap --
    def run_tleap(
        self,
        tleap_in_str: str,
        # cmd based
        if_ignore_start_up: bool = True,
        additional_search_path: List[str] = None,
        tleap_out_path: Union[str, None] = None,
        # debug
        keep_in_file: bool = False,) -> None:
        """the python wrapper of running tleap
        Args:
            tleap_in_str:
                the str content of leap.in file
            if_ignore_start_up:
                if adding the '-s' flag to ignore leaprc startup file.
            additional_search_path:
                list of addition search path of leaprc files. Each path
                in the list will be applied using a '-I' flag.
            tleap_out_path:
                file path for stdout of the tleap command. the _LOGGER level
                determines if delete the file.
            keep_in_file:
                control whether delete tleap.in file

        NOTE: run_tleap API should not handle the index alignment since it do
        not carry information of the input pdb"""
        temp_path_list = []
        # init file paths (tleap_in_path, tleap_out_path)
        fs.safe_mkdir(eh_config["system.SCRATCH_DIR"])
        tleap_in_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/tleap.in")
        temp_path_list.append(eh_config["system.SCRATCH_DIR"])
        if not keep_in_file:
            temp_path_list.append(tleap_in_path)
        if tleap_out_path is None:
            tleap_out_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/tleap.out")
            temp_path_list.append(tleap_out_path)
        # write tleap.in
        with open(tleap_in_path, "w") as of:
            of.write(tleap_in_str)
        # run tleap command
        cmd_args = f"-f {tleap_in_path} > {tleap_out_path}"
        if if_ignore_start_up:
            cmd_args = f"-s {cmd_args}"
        if additional_search_path:
            for add_path in additional_search_path:
                cmd_args = f"{cmd_args} -I {add_path}"
        try:
            self.env_manager_.run_command("tleap", cmd_args)
        except CalledProcessError as e:
            if (not e.stderr.strip()) and (not e.stdout.strip()):  # empty stderr & stdout
                # find the error information in tleap.out
                new_e = self._find_tleap_error(tleap_out_path)
                # express the error
                for e_info in new_e.error_info_list:
                    _LOGGER.error(e_info)
                raise new_e from e

        # clean up temp file if success
        fs.clean_temp_file_n_dir(temp_path_list)

    @staticmethod
    def _find_tleap_error(tleap_out_path: str) -> tLEaPError:
        """an internal used function that find the text describing the error
        from a tleap output file
        Return a tLEaPError containing all the error information"""
        error_info_pattern = r"(.*(?:FATAL|Fatal).*\n(?:.*\n)+?)(?:(?:Exiting LEaP:)|/)"
        with open(tleap_out_path) as f:
            error_info_list = re.findall(error_info_pattern, f.read())
            error_info_list = [i.strip() for i in error_info_list]
            if not error_info_list:
                _LOGGER.warning("did not found any error in tleap out file. Here is the complete output:")
                _LOGGER.warning(f.read())

        return tLEaPError(error_info_list)

    # -- antechamber/parmchk2 --
    def run_antechamber(self, in_file: str, out_file: str, net_charge: int, spin: int,
                        charge_method: str, atom_type: str = "gaff", 
                        charge_file: str=None,
                        res_name: str=None,) -> None:
        """the python wrapper of running antechamber
        Args:
            in_file: input file path
            out_file: the output molecule description file path (e.g.: mol2, prepin, ac)
            atom_type: the target atom types for generation.
                = gaff : the default
                = gaff2: for GAFF, version 2
                = amber: for PARM94/99/99SB
                = bcc : for AM1-BCC
                = sybyl: for atom types used in sybyl
            charge_method: the method for determining the atomic charge of atoms from the input file
            charge_file: the charge file. (used when charge_method="rc")
            net_charge: net_charge of the molecule. Always explicitly specify it.
            spin: multiplicity (2S+1). Always explicitly specify it.
            res_name: residue name (only used if not available in the input file)

        Unsupported options: (Add when need)
            -a additional file name
            -fa additional file format
            -ao additional file operation
                crd : only read in coordinate
                crg: only read in charge
                radius: only read in radius
                name : only read in atom name
                type : only read in atom type
                bond : only read in bond type
            -rf residue topology file name in prep input file, default is molecule.res
            -ch check file name in gaussian input file, default is molecule
            -ek QM program (mopac or sqm) keyword (in quotes); overwrites previous keywords.
            -gk gaussian keyword in a pair of quotation marks
            -gm gaussian assign memory, inside a pair of quotes, such as "%mem=1000MB"
            -gn gaussian assign number of processor, inside a pair of quotes, such as "%nproc=8"
            -gv add keyword to generate gesp file (for Gaussian 09 only) 1: yes; 0: no, the default
            -ge gaussian esp file generated by iop(6/50=1), default is g09.gesp
            -df use divcon flag, 0 - use mopac; 2 - use sqm (the default)
            -du check atom name duplications, can be yes(y) or no(n), default is yes
            -bk 4-character component Id, for ccif
            -j atom type and bond type prediction index, default is 4
                0 : no assignment
                1 : atom type
                2 : full bond types
                3 : part bond types
                4 : atom and full bond type
                5 : atom and part bond type
            -eq equalize atomic charge, default is 1 for '-c resp' and '-c bcc'
                0 : no equalization
                1 : by atomic paths
                2 : by atomic paths and geometry, such as E/Z configurations
            -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)
            -pf remove the intermediate files: can be yes (y) and no (n, default)
            -pl maximum path length to determin equivalence of atomic charges for resp and bcc.
                The smaller the value, the faster the algorithm, default is -1 (use full length),
                set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)
                -dr acdoctor mode: validate the input file a la acdoctor, yes(y, default) or no(n)"""
        in_format = self.get_file_format(in_file)
        cmd_args = ["-i", in_file,
                    "-fi", in_format,
                    "-o", out_file,
                    "-fo", self.get_file_format(out_file),
                    "-nc", str(net_charge),
                    "-m", str(spin),
                    "-at", atom_type.lower()]
        if res_name:
            cmd_args.extend(["-rn", res_name])
        # charge
        if charge_method in self.SUPPORTED_CHARGE_METHOD_MAPPER:
            charge_method = self.SUPPORTED_CHARGE_METHOD_MAPPER[charge_method]
            cmd_args.extend(["-c", charge_method])
            if charge_method == "rc":
                if charge_file:
                    cmd_args.extend(["-cf", charge_file])
                else:
                    _LOGGER.error("charge_method='rc' is used but no charge_file is provided.")
                    raise ValueError
            if charge_method == "resp":
                if in_format not in ["gesp", "gout"]:
                    _LOGGER.error(f"charge_method='resp' is used but input file format ({in_file}) is not gesp or gout.")
                    raise ValueError
        else:
            _LOGGER.error(f"found unsupported charge method {charge_method}."
                          f"(Supported keywords: {self.SUPPORTED_CHARGE_METHOD_MAPPER.keys()})"
                          "Contact author for support if it is necessary for your work.")
            raise ValueError

        self.env_manager_.run_command("antechamber", cmd_args)

        # clean up
        fs.clean_temp_file_n_dir([
            "ATOMTYPE.INF",
            "NEWPDB.PDB",
            "PREP.INF",
            "sqm.pdb",
            "sqm.in",
            "sqm.out",
        ] + glob.glob("ANTECHAMBER*"))

    def run_parmchk2(self, in_file: str, out_file: str, gaff_type: str,
                     custom_force_field: str=None,
                     is_custom_ff_type_amber: bool=True,
                     ) -> None:
        """the python wrapper of running parmchk2
        Args:
            in_file: input file path
            out_file: the output frcmod file path
            gaff_type: ff parm set, it is suppressed by "custom_force_field" option
            custom_ff: the path of a customize ff file for search.
                       (e.g.: this allows you to use ff14SB for maa)
            is_custom_ff_type_amber: is the type of custom_ff amber
            * TODO: support -c -a -w when needed"""
        cmd_args = ["-i", in_file,
                    "-f", self.get_file_format(in_file),
                    "-o", out_file,
                    "-s", gaff_type.lower(),]
        if custom_force_field:
            cmd_args.extend(["-p", custom_force_field])
            if not is_custom_ff_type_amber:
                cmd_args.extend(["-pf", "2"])

        self.env_manager_.run_command("parmchk2", cmd_args)

    # -- add_pdb --
    def run_add_pdb(self, in_prmtop: str, out_path: str, ref_pdb: str, guess: bool = False):
        """the python wrapper of running add_pdb
        Args:
            in_prmtop:
                the path of the input prmtop
            out_prmtop:
                the target output path of the modified prmtop
            ref_pdb:
                the reference pdb to add to prmtop
            guess:
                Guess atomic elements when absent from the PDB file.
                (default assumes proper element-aligned names)"""
        # run tleap command
        cmd_args = f"-i {in_prmtop} -p {ref_pdb} -o {out_path}"
        if guess:
            cmd_args = f"{cmd_args} -guess"
        self.env_manager_.run_command("add_pdb", cmd_args)

    def clean_up_add_pdb_info(self, in_prmtop: str, out_path: str):
        """remove add_pdb info in prmtop"""
        with open(in_prmtop, "r") as f, open(out_path, "w") as of:
            for line in f:
                if "RESIDUE_CHAINID" in line: # TODO move this function to PrmtopIO and make a better solution
                    break
                of.write(line)
        
    # -- cpptraj --
    def run_cpptraj(
            self,
            instr_str: str,
            log_path: Union[str, None] = None,
        ):
        """the python wrapper of running cpptraj.
        Based on https://amberhub.chpc.utah.edu/command-line-syntax/.
        Args:
            instr_str:
                the str content of cpptraj.in instruction file
            log_path:
                path for the logging from cpptraj.
                (default: {eh_config['system.SCRATCH_DIR']}/cpptraj.out)
            TODO add more when needed"""
        temp_path_list = []
        # init file paths
        fs.safe_mkdir(eh_config["system.SCRATCH_DIR"])
        in_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/cpptraj.in")
        temp_path_list.extend([eh_config["system.SCRATCH_DIR"], in_path])
        if log_path is None:
            log_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/cpptraj.out")
            temp_path_list.append(log_path)

        # write cpptraj.in
        with open(in_path, "w") as of:
            of.write(instr_str)

        # run cpptraj command
        cmd_args = f"-i {in_path} > {log_path} 2>&1"
        self.env_manager_.run_command("cpptraj", cmd_args)

        # clean up temp file if success
        fs.clean_temp_file_n_dir(temp_path_list)

    def get_n_frames_cpptraj_log(self, cpptraj_log: str) -> int:
        """get the num of frames from a cpptraj log. (this will always be in the output)"""
        with open(cpptraj_log, "r") as f:
            content = f.read()
        pattern = r"INPUT TRAJECTORIES.*\n.*\(reading 1 of (\d+)\)"
        result = int(re.search(pattern, content).group(1))
        return result

    # -- MMPBSA.py --
    def run_mmpbsa(
        self,
        # IO data
        dr_prmtop: str, dl_prmtop: str, dc_prmtop: str, sc_prmtop: str, traj_file: str,
        out_path: str, 
        # choose algorithm
        solvent_model: str,
        # algorithm config
        ion_strength: float = 0.15,
        igb: int = 5,
        fillratio: float = None, # default: 4.0
        exdi: float = None, # default: 80.0
        indi: float = None, # default: 1.0
        # operation config
        startframe: int = 1,
        endframe: int = None, # default: last frame
        interval: int = 1,
        verbose: int = 1,
        keep_files: int = 0, # NOTE seems the keep files option no longer works in terms of removing the intermediate files when prefix is used
        use_sander: bool = True,
        cluster_job_config: ClusterJobConfig = None,
        job_check_period: int = 30,
        non_armer_cpu_num: int = None,
        keep_in_file: bool = False,
        in_file: str = None,
        work_dir: str = None,
        result_in_each_frames: bool = True,
        ) -> None:
        """the python wrapper for running MMPBSA.py.MPI (MMPBSA.py not supported yet)
        Handle dispatch of using ClusterJob or not in this function.
        Args:
            dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop, traj_file:
                The input prmtop files of dry receptor, dry ligand, dry complex, solvated
                complex, and the .nc/.mdcrd file of the trajectory
            out_path:
                the path for output the .dat file
            solvent_model:
                choose between PBSA and GBSA by keywords "pbsa" or "gbsa".
            ion_strength:
                this is istrng for PB and saltcon for GB.
            startframe, endframe, interval, keep_files, igb, fillratio, verbose, 
            use_sander, exdi, indi:
                These arguments means the same as the original keywords in the .in file of
                MMPB/GBSA. See the Amber Manual (Section 34.3.1 in Amber20's case) for details.
                Only commonly used ones are included. Contact developer if you need more.
            cluster_job_config, job_check_period, non_armer_cpu_num:
                arguments for ARMer setup. 
                (see the docstring of `enzy_htp.analysis.binding.binding_energy`)
            
            keep_in_file:
                control whether clean up the .in file after finish. (force keep if in_file
                is provided by user through {in_file})
            in_file:
                this argument allows user to provide their own .in files. This will overwrite
                all the .in file related arguments.
            """
        if work_dir is None:
            work_dir = eh_config.system.SCRATCH_DIR
        fs.safe_mkdir(work_dir)
        clean_up_targets = []

        # make .in file if not provided
        if not in_file:
            in_file = fs.get_valid_temp_name(f"{work_dir}/mmpbsa.in")
            self._create_mmpbgbsa_in_file(
                in_file = in_file,
                solvent_model = solvent_model,
                ion_strength = ion_strength,
                igb = igb,
                startframe = startframe,
                endframe = endframe,
                interval = interval,
                verbose = verbose,
                keep_files = keep_files,
                use_sander = use_sander,
                fillratio = fillratio,
                exdi = exdi,
                indi = indi,
            )
        else:
            keep_in_file = True # force keep in file if it is provided.

        if cluster_job_config is not None:
            num_cores = int(cluster_job_config.node_cores)
        else:
            num_cores = non_armer_cpu_num
        # MMPBSA need num_cores <= num of frames
        num_of_frames = self.count_num_of_frames_traj(sc_prmtop, traj_file)
        if num_cores > num_of_frames:
            num_cores = num_of_frames

        if result_in_each_frames:
            temp_out_path = fs.get_valid_temp_name(f"{work_dir}/temp_mmpbsa_out.dat")
            clean_up_targets.append(temp_out_path)
        else:
            temp_out_path = None

        int_file_prefix = f"{work_dir}/_MMPBSA_"
        mmpbgbsa_command = self._make_mmpbgbsa_command(
            num_cores = num_cores,
            dr_prmtop = dr_prmtop,
            dl_prmtop = dl_prmtop,
            dc_prmtop = dc_prmtop,
            sc_prmtop = sc_prmtop,
            traj_file = traj_file,
            in_path = in_file,
            out_path = out_path,
            prefix = int_file_prefix, # control the intermediate file 
            result_in_each_frames=result_in_each_frames,
            temp_out_path=temp_out_path,
        )

        # dispatch for ARMer
        if cluster_job_config is None:
            # running locally
            mmpbgbsa_exec, mmpbgbsa_args = mmpbgbsa_command.split(" ", 1)
            this_run = self.env_manager_.run_command(
                exe=mmpbgbsa_exec,
                args=mmpbgbsa_args,
            )
            self.check_mmpbgbsa_error(out_path, this_run)
        else:
            # make job
            cluster = cluster_job_config.cluster
            res_keywords = cluster_job_config.res_keywords
            env_settings = cluster.AMBER_ENV["CPU"]
            sub_script_path = fs.get_valid_temp_name(f"{work_dir}/submit_mmpbgbsa.cmd")
            
            job = ClusterJob.config_job(
                commands=mmpbgbsa_command,
                cluster=cluster,
                env_settings=env_settings,
                res_keywords=res_keywords,
                sub_dir="./",  # because path are relative
                sub_script_path=sub_script_path)
            job.submit()
            job.wait_to_end(period=job_check_period)
            self.check_mmpbgbsa_error(out_path, job)
            clean_up_targets.extend([job.job_cluster_log, sub_script_path])
        
        # clean up
        if not keep_in_file:
            clean_up_targets.append(in_file)
        clean_up_targets.extend(glob.glob(f"{int_file_prefix}*"))
        fs.clean_temp_file_n_dir(clean_up_targets)

    def _create_mmpbgbsa_in_file(self, in_file,
            solvent_model,
            ion_strength,
            igb,
            startframe,
            endframe,
            interval,
            verbose,
            keep_files,
            use_sander,
            fillratio,
            exdi,
            indi,
        ):
        """create a mmpbgbsa .in file in {in_file}"""
        general_optional_namelist = [
            "startframe", "endframe", "interval", "verbose", "keep_files", 
            "use_sander",]
        pb_optional_namelist = [
            "fillratio", "exdi", "indi",]
        gb_optional_namelist = []

        # create in file
        in_file_lines = [
            "config generated by EnzyHTP",
        ]
        general_lines = [
            "&general",
        ]
        # general namelist
        arguments = locals()
        for k in general_optional_namelist:
            v = arguments[k]
            if v is not None:
                if isinstance(v, bool):
                    v = int(v)
                general_lines.append(f"  {k}={v},")
        general_lines.append("/")
        in_file_lines.extend(general_lines)

        if solvent_model == "pbsa":
        # pb namelist
            pb_lines = [
                "&pb",
                f"  istrng={ion_strength},",
            ]
            for k in pb_optional_namelist:
                v = arguments[k]
                if v is not None:
                    if isinstance(v, bool):
                        v = int(v)
                    pb_lines.append(f"  {k}={v},")
            pb_lines.append("/")
            in_file_lines.extend(pb_lines)

        elif solvent_model == "gbsa":
        # gb namelist
            gb_lines = [
                "&gb",
                f"  saltcon={ion_strength},",
                f"  igb={igb},",
            ]
            for k in gb_optional_namelist:
                v = arguments[k]
                if v is not None:
                    if isinstance(v, bool):
                        v = int(v)
                    gb_lines.append(f"  {k}={v},")
            gb_lines.append("/")
            in_file_lines.extend(gb_lines)
        
        else:
            _LOGGER.error(
                f"solvent_model ({solvent_model}) not supported. Supported: [pbsa, gbsa]")
            raise ValueError

        fs.write_lines(in_file, in_file_lines)

    def check_mmpbgbsa_error(
            self, 
            out_file: str, 
            stdstream_source: Union[ClusterJob, CompletedProcess, SubprocessError]
        ):
        """avoid slient error from mmpbsa run"""
        if not Path(out_file).exists():
            # collect error info
            error_info_list = []
            # 1. stdout stderr
            # error types: Fall of bus, periodic box has changed too much, illegel mem
            if isinstance(stdstream_source, ClusterJob):
                with open(stdstream_source.job_cluster_log) as f:
                    stderr_stdout = f.read()
                    error_info_list.append(f"stdout/stderr(from job log):{os.linesep*2}{stderr_stdout}")
            elif isinstance(stdstream_source, CompletedProcess):
                error_info_list.append(f"stdout:{os.linesep*2}{stdstream_source.stdout}")
                error_info_list.append(f"stderr:{os.linesep*2}{stdstream_source.stderr}")
            elif isinstance(stdstream_source, SubprocessError):
                error_info_list.append(f"stdout:{os.linesep*2}{stdstream_source.stdout}")
                error_info_list.append(f"stderr:{os.linesep*2}{stdstream_source.stderr}")
            else:
                _LOGGER.error("Only allow ClusterJob, CompletedProcess, SubprocessError as input types for `stdstream_source`")
                raise TypeError
            # 2. mdout file
            # have no experience that error pops here yet. Add when observed
            # 3. traj analysis 
            # TODO
            # Note that if_report = 0 will make traj_path = None
            # TODO give suggestion for each detected types of errors.

            _LOGGER.error(f"MMPBSA didn't finish normally.{os.linesep}{os.linesep.join(error_info_list)}")
            raise AmberMDError(error_info_list)

    def _make_mmpbgbsa_command(
            self, 
            num_cores: int,
            dr_prmtop: str, dl_prmtop: str, dc_prmtop: str, sc_prmtop: str, traj_file: str,
            in_path: str, out_path: str, temp_out_path: str, prefix: str,
            result_in_each_frames: bool = True,
        ) -> str:
        """make the mmpb/gbsa command. Overwrite existing files by default by "-O"."""
        mpi_exec = eh_config.system.get_mpi_executable(num_cores)
        mmpbsa_exec = self.config()["HARDCODE_MMMPBSA_MPI_ENGINE"]
        result = (
            f"{mpi_exec} {mmpbsa_exec} "
            f"-O "
            f"-i {in_path} "
            f"-sp {sc_prmtop} "
            f"-cp {dc_prmtop} "
            f"-rp {dr_prmtop} "
            f"-lp {dl_prmtop} "
            f"-y {traj_file} "
            f"-prefix {prefix} "
        )
        if result_in_each_frames:
            result += f"-eo {out_path} -o {temp_out_path}"
        else:
            result += f"-o {out_path} "

        return result

    def run_ante_mmpbsa(
        self,
        complex_prmtop_in: str,
        dry_complex_out: str,
        dry_receptor_out: str,
        dry_ligand_out: str,
        strip_mask: str,
        radii: str,
        receptor_mask: str = None,
        ligand_mask: str = None,
        ):
        """the python wrapper for running ante-MMPBSA.py"""
        if receptor_mask is None:
            if ligand_mask is None:
                _LOGGER.error("Must specify either 'receptor_mask' or 'ligand_mask'!")
                raise ValueError
            else:
                selection_args = ["-n", ligand_mask]
        elif ligand_mask is None:
            selection_args = ["-m", receptor_mask]
        else:
            _LOGGER.error("Must specify either 'receptor_mask' OR 'ligand_mask'! (not both) ")
            raise ValueError

        cmd_args = [
            "-p", complex_prmtop_in,
            "-c", dry_complex_out,
            "-r", dry_receptor_out,
            "-l", dry_ligand_out,
            "--radii", radii,
            "-s", strip_mask,
        ] + selection_args

        self.env_manager_.run_command("ante-MMPBSA.py", cmd_args)

    def parse_mmpbsa_result(self, mmpbsa_out_file: str, by_frames: bool) -> Dict[str, pd.DataFrame]:
        """parse the resulting mmpbsa.dat file from MMPBSA to a dictionary
        Args:
            mmpbsa_out_file: 
                the data file generated by MMPBSA.py
            by_frames: 
                specify the format of the data file to be the average one (.dat) or
                the by frames one (.csv)"""
        if by_frames:
            result = self._parse_mmpbsa_result_by_frames(mmpbsa_out_file)
        else:
            result = self._parse_mmpbsa_result_average(mmpbsa_out_file)
    
        if not result:
            _LOGGER.warning(
                f"The provided mmpbsa data file ({mmpbsa_out_file}) does not contain any data!"
            )

        return result
    
    def _parse_mmpbsa_result_by_frames(self, mmpbsa_out_file: str) -> Dict[str, pd.DataFrame]:
        """parse the result when -eo is specified for the output"""
        gb_pattern = (
            r"GENERALIZED BORN:(?:.|\n)+?DELTA Energy Terms\n"
            r"((?:.|\n)+?)"
            r"\n\n"
        )
        pb_pattern = (
            r"POISSON BOLTZMANN:(?:.|\n)+?DELTA Energy Terms\n"
            r"((?:.|\n)+?)"
            r"\n\n"
        )
        result = {}

        with open(mmpbsa_out_file) as f:
            f_str = f.read()
        # GB
        gb_search = re.search(gb_pattern, f_str)
        if gb_search is not None:
            gb_table = gb_search.group(1).strip()
            gb_result = pd.read_csv(StringIO(gb_table))
            result["gbsa"] = gb_result
        # PB
        pb_search = re.search(pb_pattern, f_str)
        if pb_search is not None:
            pb_table = pb_search.group(1).strip()
            pb_result = pd.read_csv(StringIO(pb_table))
            result["pbsa"] = pb_result

        return result

    def _parse_mmpbsa_result_average(self, mmpbsa_out_file: str) -> Dict[str, pd.DataFrame]:
        gb_pattern = (
            r"GENERALIZED BORN:(?:.|\n)+?Differences \(Complex - Receptor - Ligand\):\n"
            r"((?:.|\n)+?)"
            r"-------------------------------------------------------------------------------\n"
            r"-------------------------------------------------------------------------------"
        )
        pb_pattern = (
            r"POISSON BOLTZMANN:(?:.|\n)+?Differences \(Complex - Receptor - Ligand\):\n"
            r"((?:.|\n)+?)"
            r"-------------------------------------------------------------------------------\n"
            r"-------------------------------------------------------------------------------"
        )
        result = {}

        with open(mmpbsa_out_file) as f:
            f_str = f.read()
        # GB
        gb_search = re.search(gb_pattern, f_str)
        if gb_search is not None:
            gb_table = gb_search.group(1).strip()
            gb_result = self._parse_pb_gb_table(gb_table)
            result["gbsa"] = gb_result
        # PB
        pb_search = re.search(pb_pattern, f_str)
        if pb_search is not None:
            pb_table = pb_search.group(1).strip()
            pb_result = self._parse_pb_gb_table(pb_table)
            result["pbsa"] = pb_result

        return result
    
    def _parse_pb_gb_table(self, table_str: str) -> pd.DataFrame:
        """the table looks like this:
        Energy Component            Average              Std. Dev.   Std. Err. of Mean
        -------------------------------------------------------------------------------
        VDWAALS                    -20.3834                3.1308              0.3131
        EEL                        -19.1765                7.2973              0.7297
        EGB                         20.8964                4.2728              0.4273
        ESURF                       -3.8046                0.1256              0.0126

        DELTA G gas                -39.5599                8.4658              0.8466
        DELTA G solv                17.0918                4.2257              0.4226

        DELTA TOTAL                -22.4680                5.2186              0.5219
        Returns: (like this) pd.DataFrame
                          mean      sd     sem
            VDWAALS      -20.3834  3.1308  0.3131
            EEL          -19.1765  7.2973  0.7297
            EGB           20.8964  4.2728  0.4273
            ESURF         -3.8046  0.1256  0.0126
            DELTA G gas  -39.5599  8.4658  0.8466
            DELTA G solv  17.0918  4.2257  0.4226
            DELTA TOTAL  -22.4680  5.2186  0.5219"""

        data_pattern = r"([A-z]+(?: [A-z]+)*)( *[0-9\-\.]+)( *[0-9\-\.]+)( *[0-9\-\.]+)"
        data_lines = re.findall(data_pattern, table_str)
        data_dict = {}
        for line in data_lines:
            data_dict[line[0].strip()] = [float(line[1].strip()), float(line[2].strip()), float(line[3].strip())]
        
        result_df = pd.DataFrame.from_dict(
            data_dict, 
            orient="index",
            columns=["mean", "sd", "sem"]
            )
        
        return result_df

    # -- parmed --
    def run_parmed(
        self,
        parmed_in_str: str,
        prmtop_file: str,
        inpcrd_file: str = None,
        log_file: str = None,
        overwrite: bool = True,
        ):
        """the python wrapper for running parmed"""
        temp_dir = eh_config.system.SCRATCH_DIR
        fs.safe_mkdir(temp_dir)
        temp_in_file = fs.get_valid_temp_name(f"{temp_dir}/temp_parmed.in")
        with open(temp_in_file, "w") as of:
            of.write(parmed_in_str)

        cmd_args = [
            "-i", temp_in_file,
            "-p", prmtop_file,
        ]
        if inpcrd_file is not None:
            cmd_args += ["-c", inpcrd_file]
        if log_file is not None:
            cmd_args += ["-l", log_file]
        if overwrite:
            cmd_args += ["-O"]

        self.env_manager_.run_command("parmed", cmd_args)

        # clean up
        fs.clean_temp_file_n_dir([
            temp_in_file
        ])

    # -- index mapping --
    def get_amber_index_mapper(self, stru: Structure) -> Dict[str, Dict[Union[Residue, Atom], Union[int, tuple]]]:
        """get a mapper for objects in {stru} and in Amberilzed PDB of stru
        Use tLeap to get a PDB and align indexes.
        Returns:
            {
            "residue" : {Residue(): amber_1_index_residue_idx_1, ...},
            "atom" : {Atom(): amber_1_index_atom_idx_1}
            }

        This will work as long as the order of residues remains the same after tleap process and atom
        have unique names in each residue.
        TODO solve this for the 1Q4T case."""
        # init files
        temp_dir = eh_config['system.SCRATCH_DIR']
        fs.safe_mkdir(temp_dir)
        temp_pdb = fs.get_valid_temp_name(f"{temp_dir}/amber_index_mapping_source.pdb")
        temp_amber_pdb = fs.get_valid_temp_name(f"{temp_dir}/amber_index_mapping_amber.pdb")
        # init stru
        pdb_io.PDBParser().save_structure(
            temp_pdb, stru,
            if_renumber=False, if_fix_atomname=False)
        self.tleap_clean_up_stru(
            input_pdb=temp_pdb,
            out_path=temp_amber_pdb,
            if_align_index=False,
        )
        # deduce mapping
        ref_stru = pdb_io.PDBParser().get_structure(temp_amber_pdb)
        ref_mapper = stru.get_residue_index_mapper(ref_stru)
        result = {
            "residue" : {},
            "atom" : {},
        }
        ref_reskey_mapper = ref_stru.residue_mapper
        for k, res in stru.residue_mapper.items():
            amber_key = ref_mapper[k]
            result["residue"][res] = amber_key
            for atom in res.atoms:
                amber_atom_idx = ref_reskey_mapper[amber_key].find_atom_name(atom.name).idx
                result["atom"][atom] = amber_atom_idx
        fs.clean_temp_file_n_dir([
            temp_dir,
            temp_pdb,
            temp_amber_pdb,
        ])
        # TODO add a cache mechanism
        return result

    def get_amber_atom_index(self, atoms: List[Atom]) -> List[int]:
        """get atom index in an Amberized PDB of the Structre() containing these Atom()s"""
        stru = atoms[0].root()
        aid_mapper = self.get_amber_index_mapper(stru)
        return [aid_mapper["atom"][at] for at in atoms]

    def rename_atoms(self, stru: Structure) -> None: # TODO(high piror) move to structure_io https://github.com/ChemBioHTP/EnzyHTP/pull/162#discussion_r1473217587
        """Renames residues and atoms to be compatible with Amber naming and functions.
        
        Args:
            stru: The Structure() to perform renaming on.

        Returns:
            Nothing. 
        """
        nterm_mapper:Dict = {"H1":"1H", "H2":"2H", "H3":"3H"}
        leu_mapper:Dict = {
            "1HD1":"HD11",
            "2HD1":"HD12",
            "3HD1":"HD13",
            "1HD2":"HD21",
            "2HD2":"HD22",
            "3HD2":"HD23",
            }
        _LOGGER.info("Beginning renaming...")
        changed_residues:int = 0
        changed_atoms:int = 0
        for res in stru.residues:
            if not res.is_canonical():
                continue
            

            if res.name == 'LEU':
                for aa in res.atoms:
                    if aa.name in leu_mapper:
                        aa.name = leu_mapper[aa.name]
                        changed_atoms += 1                            
        
        _LOGGER.info(f"Finished renaming! Changed {changed_residues} residues and {changed_atoms} atoms.")

    def get_amber_residue_index(self, residues: List[Atom]) -> List[int]:
        """"""
        raise Exception("TODO")

    def get_amber_chain_index(self, chains: List[Atom]) -> List[str]:
        """"""
        raise Exception("TODO")

    # -- amber mask --
    # region Amber Mask I/O
    def get_amber_mask(self, target: StruSelection, reduce: bool) -> str:
        """get the Amber mask based on a StruSelection() described by the mask
        Args:
            target:
                the target StruSelection() described by the mask
            reduce:
                whether the function reduces the mask for better readability
        Ruturns:
            the Amber Mask (e.g.: "'@C,CA,N'")
        Notes:
            reference: https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/
            There are several ways to represent the same set of Atom()s. If
            not reduced, we use directly the atom numlist. It will be Structure()
            specific and do not respond to mutation as Structure() is not the same."""
        # san check
        target.check_consistent_topology()

        atom_idxes = self.get_amber_atom_index(target.atoms)
        reduced_idxes = get_interval_str_from_list(atom_idxes)
        mask = f"'@{reduced_idxes}'"

        if reduce:
            pass 
            # TODO. reduce based on residue affilation etc. 
            # use methods from StruSelection

        return mask
            
    # endregion

    # -- sander/pmemd --
    # region - mdin/mdout file I/O -
    def read_from_mdin(self, mdin_file_path: str) -> Dict:
        """the knowledge of the mdin format as the input of sander/pmemd.
        return a dictionary that use standard keys according to AmberMDStep.md_config_dict"""
        _LOGGER.error("dont support this yet. contact developer if you have to use this.")
        # TODO make raw_dict a class when work on this.
        raise Exception
        # TODO
        # 1. parse to dict
        _read_from_mdin_to_raw_dict()
        # 2. mapper keywords and build objects

    def _read_from_mdin_to_raw_dict():
        """read from mdin file and parse it to a raw dictionary
        where key is Amber keywords as is. (defined in
        _write_to_mdin_from_raw_dict())"""

    def write_disang_file(self, raw_rs_dict_list: List[dict], rs_path: str):
        """write {raw_rs_dict} to {rs_path}"""
        disang_lines = []
        for raw_dict in raw_rs_dict_list:
            cons_lines = [
                "&rst"
            ]
            for k, v in raw_dict.items():
                if k == "iat" or k.startswith("igr"):
                    v = map(str, v)
                    cons_lines.append(f" {k}={','.join(v)},")
                else:
                    cons_lines.append(f" {k}={v},")
            cons_lines = cons_lines + [
                "&end",
            ]
            disang_lines.extend(cons_lines)
        disang_lines.append("")
        fs.write_lines(rs_path, disang_lines)

    REDIRECTION_FILE_WRITER = {
        "DISANG" : write_disang_file,
    }

    def write_to_mdin(self, md_config_dict: Dict, out_path: str, path_rel_to: str = None):
        """parse a AmberMDStep.md_config_dict to a mdin file."""
        # parse data_dict
        raw_data_dict = self._parse_md_config_dict_to_raw(md_config_dict)
        # write to file
        self._write_to_mdin_from_raw_dict(raw_data_dict, out_path, path_rel_to=path_rel_to)

    def _parse_md_config_dict_to_raw(self, md_config_dict: Dict) -> Dict:
        """parse AmberMDStep.md_config_dict to the raw data dict for mdin
        hypothesis: if file_redirection is involved a wt type=END is needed?

        Return: the raw data dict.
        TODO should make the raw dict a class to strictly define it."""
        wt_list = []
        file_redirection_dict = {}
        group_info_list = []

        # imin, ntx, irest, ntc, ntf
        imin = int(md_config_dict["minimize"])
        ntx, irest = self.MD_RESTART_MAPPER[md_config_dict["restart"]]
        ntc, ntf = self.MD_TIMESTEP_SHAKE_MAPPER[md_config_dict["timestep"]]

        # region == imin variant == (some keyword dont work when imin=1 or 0 each)
        imin_or_not_cntrl = {
            "imin" : imin,
        }

        if imin == 0: # MD
            # nstlim
            timestep = md_config_dict["timestep"]
            dt = timestep * 1000
            raw_nstlim = md_config_dict["length"] / timestep
            nstlim = mh.round_by(raw_nstlim, 0.5)
            if not mh.is_integer(raw_nstlim):
                _LOGGER.warning(f"length ({md_config_dict['length']}) is not divisible by"
                                f" timestep ({md_config_dict['timestep']}). Setting nstlim"
                                f" to {nstlim}. (i.e.:actual length={nstlim * timestep} ns)")

            # ntpr
            ntpr = max(int(self.config()["HARDCODE_NTPR_RATIO"] * nstlim), 1)

            # ntwx
            ntwx = 0
            if md_config_dict["if_report"]:
                ntwx = mh.round_by(md_config_dict["record_period"] / timestep, 0.5)

            # region == temperture related keywords ==
            temperature_data = md_config_dict["temperature"]
            if isinstance(temperature_data, float):
                temp_cntrl = {'temp0': temperature_data}
            elif isinstance(temperature_data, list):
                temp_wt_list = []
                temp_cntrl = {'tempi': temperature_data[0][1], 'temp0': temperature_data[-1][1]}
                temperature_data = [ (mh.round_by(x/timestep, 0.5), y) for x,y in temperature_data]
                for (step0, temp0), (step1, temp1) in mh.get_section_from_endpoint(temperature_data):
                    temp_wt_list.append({
                        'type': 'wt',
                        'config': {
                            'type': "'TEMP0'",
                            'istep1': step0, 'istep2': step1,
                            'value1': temp0, 'value2': temp1,}
                        })
                wt_list.extend(temp_wt_list)
            else:
                _LOGGER.error("temperature can only be a float or a list")
                raise TypeError
            # endregion

            # ntt
            ntt = self.MD_THERMOSTAT_MAPPER.get(md_config_dict["thermostat"], None)
            if ntt is None:
                _LOGGER.error(f"thermostat keyword: {md_config_dict['thermostat']} is not supported."
                            f" Supported list: {self.MD_THERMOSTAT_MAPPER.keys()}")
                raise ValueError
            ntt_cntrl = {'ntt': ntt}
            if ntt == 3:
                ntt_cntrl.update({'gamma_ln': self.config()["HARDCODE_GAMMA_LN"],})
            
            # ntp
            if imin == 1:
                ntp_cntrl = {}
            else:
                ntp = self.MD_PRESSURE_SCALING_MAPPER.get(md_config_dict["pressure_scaling"], None)
                if ntp is None:
                    _LOGGER.error(f"pressure_scaling keyword: {md_config_dict['pressure_scaling']} is not supported."
                                f" Supported list: {self.MD_PRESSURE_SCALING_MAPPER.keys()}")
                    raise ValueError
                ntp_cntrl = {"ntp" : ntp}
            
            # ntb
            if imin == 1:
                ntb_cntrl = {}
            elif ntp == 0:
                ntb_cntrl = {"ntb" : 1}
            else:
                ntb_cntrl = {"ntb" : 2}

            imin_or_not_cntrl = imin_or_not_cntrl | {
                'nstlim': nstlim, 'dt': dt,
                } | temp_cntrl | ntt_cntrl | ntb_cntrl | ntp_cntrl | {
                'iwrap': self.config()["HARDCODE_IWRAP"],
                'ig': self.config()["HARDCODE_IG"],
                }

        else: # minimization
            # maxcyc, ncyc
            maxcyc = md_config_dict["length"]
            ncyc = max(mh.round_by(maxcyc * self.config()["HARDCODE_NCYC_RATIO"], 0.5), 1)
            # ntpr
            ntpr = max(int(self.config()["HARDCODE_NTPR_RATIO"] * maxcyc), 1)
            # ntwx
            ntwx = 0

            imin_or_not_cntrl = imin_or_not_cntrl | {
                "maxcyc" : maxcyc,
                "ncyc" : ncyc,
            }
        # endregion

        # region == constraint ==
        constraints = md_config_dict["constrain"]
        cart_freeze = []
        geom_cons = []
        if constraints:
            for cons in constraints:
                cons: StructureConstraint
                if cons.constraint_type in self.config()["SUPPORTED_CONSTRAINT_TYPE"]:
                    if cons.is_cartesian_freeze():
                        cart_freeze.append(cons)
                    else:
                        geom_cons.append(cons)
                else:
                    _LOGGER.error(
                        f"Dont support {cons.constraint_type} yet in AmberInterface. Supported: {self.config()['SUPPORTED_CONSTRAINT_TYPE']}")
                    raise ValueError

        # ntr and cartesian restraint
        ntr_cntrl = {}
        if cart_freeze:
            cons = merge_cartesian_freeze(cart_freeze)
            restraintmask = self.get_restraintmask(cons, reduce=True)
            restraint_wt = cons.params["amber"]["restraint_wt"]
            ntr_cntrl = {"ntr": 1,
                         'restraint_wt': restraint_wt,
                         'restraintmask': restraintmask,}

        # nmropt
        nmropt_cntrl = {}
        if geom_cons:
            geom_cons: List[StructureConstraint]
            nmropt_cntrl = {'nmropt': 1}
            # figure out path
            disang_path: str = geom_cons[0].params["amber"]["rs_filepath"]
            if disang_path.startswith("{mdstep_dir}"):
                disang_path = disang_path.lstrip("{mdstep_dir}")
                disang_path = f"{md_config_dict['mdstep_dir']}/{disang_path}"
            disang_path = fs.get_valid_temp_name(disang_path)
            # figure out content
            disang_content_list = []
            for cons in geom_cons:
                raw_rs_dict = self._parse_cons_to_raw_rs_dict(cons)
                disang_content_list.append(raw_rs_dict)
            # assemble
            nmropt_file_redirection = {
                "DISANG" : {
                    "path": disang_path, "content": disang_content_list}
            }
            file_redirection_dict.update(nmropt_file_redirection)
        # endregion

        # determine the &wt END
        if wt_list or file_redirection_dict:
            wt_list.append(
                {'type': 'wt',
                 'config': {'type': "'END'",}},
                )

        # assemble namelists
        namelists = [
            {'type': 'cntrl',
            'config': imin_or_not_cntrl | {
                'ntx': ntx, 'irest': irest,
                'ntc': ntc, 'ntf': ntf,
                'cut': self.config()["HARDCODE_CUT"],
                'ntpr': ntpr, 'ntwx': ntwx,
                } | ntr_cntrl | nmropt_cntrl | {
            }},
        ] + wt_list

        # compile the raw dict using md_config_dict
        raw_dict = {
            'title': md_config_dict["name"],
            'namelists': namelists,
            'file_redirection': file_redirection_dict,
            'group_info': group_info_list,
        }

        return raw_dict

    def _parse_cons_to_raw_rs_dict(self, cons: StructureConstraint) -> Dict:
        """parse StructureConstraint() to a raw dict for writing the DISANG file.
        Expression containing value are also parsed in this function.
        (based on section 27.1 in Amber20 manual)
        Example rs dict
        {
        'iat':[100,101],
        'r1': 0.50,'r2': 0.96,'r3': 2.50,'r4': 3.50,
        'rk2':0.0,'rk3':100.0,'ir6':1,'ialtd':0
        }"""
        target_keys = "r1 r2 r3 r4 rk2 rk3 ir6 ialtd ifvari".split() # TODO add more when needed
        source_dict = cons.params["amber"]
        # iat & igr
        if cons.is_group_constraint():
            atom_idxes = [-1] * cons.num_groups
            result = {"iat": atom_idxes}
            for i, grp_atom in enumerate(cons.atoms_by_groups):
                result[f"igr{i+1}"] = self.get_amber_atom_index(grp_atom)
        else:
            atom_idxes = self.get_amber_atom_index(cons.atoms)
            result = {"iat": atom_idxes}

        # other keys
        for k in target_keys:
            if k in source_dict:
                v = source_dict[k] 
    
                # parse x containing value
                if k in "r1 r2 r3 r4".split():
                    if isinstance(v, str) and v.count("x") > 0:
                        expr = sympify(v)
                        v = round(
                                float(expr.subs("x", cons.target_value)),
                                ndigits=3
                            )

                result[k] = v
        
        return result

    def get_restraintmask(self, cons: CartesianFreeze, reduce: bool= True) -> str:
        """get the str of restraintmask from CartesianFreeze().
        Args:
            cons: the target constriant
            reduce: whether the function try to reduce the format of the mask
        Returns:
            the restraintmask"""
        if reduce and cons.constraint_type == "backbone_freeze":
            return "'@C,CA,N'"
        else:
            result = self.get_amber_mask(
                StruSelection(cons.atoms),
                reduce=reduce)
            return result

    def _write_to_mdin_from_raw_dict(self, raw_dict: Dict, out_path: str, path_rel_to: str = None):
        """write to mdin file from a raw dictionary
        where key is Amber keywords.
        Based on section 19.5 of Amber20 manual.
        Format of raw_dict:
            {'title': 'heat',
            'namelists': [
                {'type': 'cntrl',
                'config': {'k': v, ...}},
                {'type': 'wt',
                'config': {'k': v, ...}},],
            'file_redirection': {
                'DISANG': './MD/0.rs'},
            'group_info': [],}
            The 4 main keys are required.
        (see example in test_amber_interface.py::test_write_to_mdin_from_raw_dict)
        
        path_rel_to: 
            if this argument is specified, the file path inside of the .in file
            will be written as the relative path to {path_rel_to}. However, the
            path stored in python variables are still absolute or relative to cwd
            so that they can be valid."""
        # title
        mdin_lines: List[str] = [
            f"{raw_dict['title']}",
        ]
        # namelists
        for namelist in raw_dict["namelists"]:
            mdin_lines.append(f" &{namelist['type']}")
            for k, v in namelist['config'].items():
                mdin_lines.append(f"  {k} = {v},")
            mdin_lines.append(" /")
        # file redirection
        for k, v in raw_dict["file_redirection"].items():
            v_path = v["path"]
            v_path = self._reduce_path_in_mdin(v_path)
            v_content = v["content"]
            # add the redirection line
            if path_rel_to:
                mdin_lines.append(f" {k}={fs.relative_path_of(v_path, path_rel_to)}")
            else:
                mdin_lines.append(f" {k}={v_path}")
            # write the redirection file
            if v_content:
                self.REDIRECTION_FILE_WRITER[k](self, v_content, v_path)
        # group information
        for group_info in raw_dict["group_info"]: #TODO cant find any example of how this looks like
            _LOGGER.error("group info support not developed yet. contact author if you have to use this")
            raise ValueError
        # final empty line
        mdin_lines.append("")

        # security check
        for line in mdin_lines:
            if len(line.strip()) > 80:
                _LOGGER.warning("found mdin line exceed 80 characters.\n"
                                f"{line.strip()}\n"
                                " PMEMD may not be happy about this."
                                " Consider this bug from Amber if you get 'STOP PMEMD Terminated Abnormally'")

        fs.write_lines(out_path, mdin_lines)

    def _reduce_path_in_mdin(self, target_path: str) -> str:
        """compress the paths that >70 characters to be shorted than 70 characters.
        (we dont use 80 because the line limitation also counts in the variable name like
        DISANG etc.)"""
        if len(target_path) <= 70:
            return target_path
        else:
            _LOGGER.info("found file path exceeding length limit of mdin. reducing it.")
            # 0. resolve path
            target_path: Path = Path(target_path)
            target_path = target_path.resolve()

            # 1. relative path
            result = target_path.relative_to(Path.cwd())
            if len(str(result)) <= 70:
                return str(result)

            # 2. symlink in SCRATCH
            _LOGGER.debug("1st attempt using relative path failed. making it symlink.")
            result = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/{target_path.name}.lnk", is_symlink=True)
            if len(result) <= 70:
                fs.safe_mkdir(eh_config['system.SCRATCH_DIR'])
                Path(result).symlink_to(target_path)
                return result
            
            # 3. symlink in SCRATCH relative
            _LOGGER.debug("2nd attempt failed. making the symlink using relative path.")             
            result = Path(result).relative_to(Path.cwd())
            if len(str(result)) <= 70:
                fs.safe_mkdir(eh_config['system.SCRATCH_DIR'])
                Path(result).symlink_to(target_path)
                return str(result)

            # 4. symlink in ./
            _LOGGER.debug("3nd attempt failed. making the symlink under cwd.")
            result = fs.get_valid_temp_name(f"./{target_path.name}.lnk", is_symlink=True)
            if len(result) <= 70:
                Path(result).symlink_to(target_path)
                return result
            
            # 5. shortest symlink that have to work
            result = fs.get_valid_temp_name(f"./symlink.lnk", is_symlink=True)
            _LOGGER.warning(f"4nd attempt failed. making {result} for {target_path}.")
            Path(result).symlink_to(target_path)
            return result

    def read_from_mdout(self, mdout_file_path: str) -> Dict:
        """the knowledge of the mdout format as the log of sander/pmemd.
        return a dictionary."""
        _LOGGER.error("dont support this yet. contact developer if you have to use this.")
        raise Exception
        # TODO
        # 1. parse to dict
        # 2. mapper keywords and build objects
    # endregion

    def get_md_executable(self, core_type: str, num_cores: int) -> str:
        """get the name of the md executable depending on core_type
        return the executable and mpi prefix if needed."""
        result = self.config()["HARDCODE_MD_ENGINE"][core_type]
        if core_type == "cpu":
            
            mpi_exec = eh_config._system.get_mpi_executable(num_cores)
            result = f"{mpi_exec} {result}"
        return result

    def is_md_completed(self, mdout_file: str) -> bool:
        """check whether an Amber MD run is completed.
        clues of incomplete:
        1. no mdout file
        2. no TIMINGS section in mdout file"""
        if os.path.exists(mdout_file):
            with open(mdout_file) as f:
                return "5.  TIMINGS" in f.read() # TODO add more logic when needed.
        else:
            return False

    # def run_sander/pmemd() >>> see AmberMDStep.run() or .make_job()

    # == engines ==
    # (engines for science APIs)
    # -- mutation/clean up --
    def tleap_clean_up_stru(
        self,
        input_pdb: str,
        out_path: str,
        if_align_index: bool = True,
        amber_lib: str = "leaprc.protein.ff14SB",) -> None:
        """Method that uses tleap to clean up {input_pdb} by loading and saving the PDB with
        the {amber_lib}.
        Typical changes are:
            - complete missing atoms based on the AA defination in {amber_lib}
            of the correponding residue name
            - change terminal atom names
            - (if_align_index = False)renumber residue index and atom index from 1 and
            ignore chain index
            - add TER around all residue unit that is not recognized or is a ligand.
            (The might be harmful when a structure have a modified animo acid TODO)

        Args:
            input_pdb:
                the path of the input pdb for cleaning
            out_path:
                the path of the output pdb after cleaning
            if_align_index:
                whether or not align index back after cleaning.
            amber_lib:
                the amber library used for cleaning.

        Application:
            Used in mutation.api.mutate_stru_with_tleap()"""

        tleap_in_lines: List[str] = [
            f"source {amber_lib}",
            f"a = loadpdb {input_pdb}",
            f"savepdb a {out_path}",
            "quit",
        ]
        tleap_in_str = "\n".join(tleap_in_lines)
        self.run_tleap(tleap_in_str)
        if if_align_index:
            # temp file
            temp_dir = eh_config["system.SCRATCH_DIR"]
            fs.safe_mkdir(temp_dir)
            renumbered_pdb = fs.get_valid_temp_name(f"{temp_dir}/tleap_out_renumbered.pdb")
            # renumber
            index_mapper = pdb_io.get_index_mapper_from_pdb(out_path, input_pdb, method="by_order")
            pdb_io.restore_pdb_index(index_mapper, out_path, renumbered_pdb)
            shutil.move(renumbered_pdb, out_path)
            fs.clean_temp_file_n_dir([temp_dir, renumbered_pdb])

    # -- MD --
    def build_md_parameterizer(
            self,
            # find default values in AmberConfig
            # method
            force_fields: List[str] = "default",
            charge_method: str = "default",
            resp_engine: str = "default",
            resp_lvl_of_theory: str = "default",
            # workflow
            ncaa_param_lib_path: str = "default",
            force_renew_ncaa_parameter: bool = "default",
            ncaa_net_charge_engine: str = "default",
            ncaa_net_charge_ph: float = "default",
            solvate_box_type: str = "default",
            solvate_box_size: float = "default",
            # mod_aa_engine TODO ask if the other way really work
            # MCPY.py related TODO add when use
            # GB related
            gb_radii: int = None,
            # temp paths
            parameterizer_temp_dir: str = "default", # default: {config[system.SCRATCH_DIR]}/amber_parameterizer
            additional_tleap_lines: List[str] = None,
            # debug options
            keep_tleap_in: Union[bool, str] = "default",
            ) -> AmberParameterizer:
        """the constructor for AmberParameterizer
        Args:
            force_fields:
                The list of force fields used for parameterization in Amber tleap format
                (e.g.: ["leaprc.protein.ff14SB", "leaprc.gaff", "leaprc.water.tip3p"])
                (NOTE that if the user has access to AmberTools24 or later. It is recommend
                to use `abcg2` charge method for gaff2 and `bcc` for gaff according to Page
                312 of Amber24 manual.)
            charge_method:
                The method used for determine the atomic charge.
                This method is applied to parameterization of ligand, modified AA,
                and metal binding site.
            (the following args only apply when charge_method="RESP")
                resp_engine:
                    The engine for calculating the RESP charge.
                resp_lvl_of_theory:
                    The level of theory for calculating the RESP charge
            ncaa_param_lib_path:
                The path of the non-CAA parameter library. This is where all generated
                NCAA params goes to. It will prevent redundant generation of same NCAAs.
                Normally we suggest setting this to a directory that contains all workflows
                of a same wild-type/template enzyme.
                * The NCAA-file correspondence is determined by the 
                  (1) 3-letter name in the file
                  (2) (if not 1 not exist) the file name
                * Setting this to a path that is too general may cause conflict when different
                NCAAs have the same name. (e.g.: different tautomer or general res name like LIG)
            force_renew_ncaa_parameter:
                Whether force renew the parameter files (frcmod etc.) for all ncaa
                (ligand, mod AA, or metal)
            ncaa_net_charge_engine:
                The engine the determines the net charge of NCAA if none is assigned in NCAA objects (Ligand,
                ModifedResidue TODO, MetalUnit)
            ncaa_net_charge_ph:
                The pH value used in determining the net charge of NCAA.
            solvate_box_type:
                The shape of the solvation box.
            solvate_box_size:
                The size of the solvation box.
            gb_radii:
                The igb number - the effective GB radii used in the Generalized Born calculation.
                This will influence the GB radii in the prmtop file and are only used implicit solvent calculations.
            parameterizer_temp_dir: (default: {config[system.SCRATCH_DIR]}/amber_parameterizer)
                The temporary working directory that contains all the files generated by the AmberParameterizer
            additional_tleap_lines:
                handle for adding additional tleap lines before generating the parameters.
            keep_tleap_in:
                whether keep tleap.in file generated during the process. (default: False)"""
        # tool: write the below code
        # print(AmberInterface._generate_default_assigning_lines_for_build_md_parameterizer_or_step(locals().items()))

        # init default values
        type_hint_sticker: AmberConfig
        if force_fields == "default":
            force_fields = self.config()["DEFAULT_FORCE_FIELDS"]        
        if charge_method == "default":
            charge_method = self.config()["DEFAULT_CHARGE_METHOD"]
        if resp_engine == "default":
            resp_engine = self.config()["DEFAULT_RESP_ENGINE"]
        if resp_lvl_of_theory == "default":
            resp_lvl_of_theory = self.config()["DEFAULT_RESP_LVL_OF_THEORY"]
        if ncaa_param_lib_path == "default":
            ncaa_param_lib_path = self.config()["DEFAULT_NCAA_PARAM_LIB_PATH"]
        if force_renew_ncaa_parameter == "default":
            force_renew_ncaa_parameter = self.config()["DEFAULT_FORCE_RENEW_NCAA_PARAMETER"]
        if ncaa_net_charge_engine == "default":
            ncaa_net_charge_engine = self.config()["DEFAULT_NCAA_NET_CHARGE_ENGINE"]
        if ncaa_net_charge_ph == "default":
            ncaa_net_charge_ph = self.config()["DEFAULT_NCAA_NET_CHARGE_PH"]
        if solvate_box_type == "default":
            solvate_box_type = self.config()["DEFAULT_SOLVATE_BOX_TYPE"]
        if solvate_box_size == "default":
            solvate_box_size = self.config()["DEFAULT_SOLVATE_BOX_SIZE"]
        if parameterizer_temp_dir == "default":
            parameterizer_temp_dir = self.config()["DEFAULT_PARAMETERIZER_TEMP_DIR"]
        if keep_tleap_in == "default":
            keep_tleap_in = self.config()["DEFAULT_KEEP_TLEAP_IN"]

        return AmberParameterizer(
            self,
            force_fields,
            charge_method,
            resp_engine,
            resp_lvl_of_theory,
            ncaa_param_lib_path,
            force_renew_ncaa_parameter,
            ncaa_net_charge_engine,
            ncaa_net_charge_ph,
            solvate_box_type,
            solvate_box_size,
            gb_radii,
            parameterizer_temp_dir,
            additional_tleap_lines,
            keep_tleap_in,)

    @staticmethod
    def _generate_default_assigning_lines_for_build_md_parameterizer_or_step(parameter_mapper: dict) -> str:
        """function only used for generate the code of the init default value section
        of build_md_parameterizer"""
        result = ""
        for k, v in parameter_mapper:
            if v == "default":
                result += (f"""if {k} == "default":
            {k} = self.config()["DEFAULT_{k.upper()}"]
        """)
        return result

    def antechamber_ncaa_to_moldesc(self,
                                    ncaa: NonCanonicalBase,
                                    out_path: Union[str, None] = None,
                                    gaff_type: str = "GAFF",
                                    charge_method: str = "AM1BCC",
                                    cluster_job_config: Dict=None,) -> str:
        """use antechamber to generate .mol2/.ac file for ligand/modified amino acid.
        Args:
            ncaa: the target Ligand/ModifiedAminoAcid
            out_path: the path of the output molecule description file. (use ext here to determine target format)
            gaff_type: the ff type used for NCAA. This influence the atom type in the moldesc file
            charge_method: the method for generation of atomic charges
            cluster_job_config: the configuration for submitting cluster jobs (used for Gaussian if charge_method='resp')
        Return:
            the out_path
            """
        # san check
        if (ncaa.multiplicity is None) or (ncaa.net_charge is None):
            _LOGGER.error(f"supplied NCAA ({ncaa.name}) does not have charge and spin."
                          " ALWAYS check and explicit assign it using"
                          " Structure.assign_ncaa_chargespin()")
            raise ValueError

        if ncaa.is_modified():
            atom_type = "AMBER"
        else:
            atom_type = gaff_type

        multiplicity = ncaa.multiplicity
        net_charge = ncaa.net_charge

        # init_path
        if out_path is None:
            if ncaa.is_modified():
                out_path = fs.get_valid_temp_name(f"{eh_config['system.NCAA_LIB_PATH']}/{ncaa.name}_{charge_method}-{atom_type}.ac")
            else:
                out_path = f"{eh_config['system.NCAA_LIB_PATH']}/{ncaa.name}_{charge_method}-{gaff_type}.mol2"

        # 1. make ligand PDB
        temp_dir = eh_config["system.SCRATCH_DIR"]
        fs.safe_mkdir(temp_dir)
        temp_pdb_path = fs.get_valid_temp_name(f"{temp_dir}/{ncaa.name}.pdb")
        if ncaa.is_modified():
            # 1.1. Capping - cap C-terminal with OH and N-terminal with H
            ncaa_region = create_region_from_residues(residues=[ncaa], nterm_cap="H", cterm_cap="OH")
            ncaa = ncaa_region.convert_to_structure(cap_as_residue=False)
        pdb_io.PDBParser().save_structure(temp_pdb_path, ncaa)
        input_file = temp_pdb_path

        # 2. determine charge
        charge_method = self.SUPPORTED_CHARGE_METHOD_MAPPER.get(charge_method, None)
        if charge_method:
            if charge_method == "bcc":
                pass
            if charge_method == "resp":
            # 2.1. Run RESP calculation (option)
                gout_file = self._calculate_charge_with(charge_method, ncaa, cluster_job_config)
                input_file = gout_file
        else:
            _LOGGER.error(f"found unsupported charge method {charge_method}."
                          f"(Supported keywords: {self.SUPPORTED_CHARGE_METHOD_MAPPER.keys()})")
            raise ValueError
        

        # 3. run antechamber on the PDB get mol2
        self.run_antechamber(in_file=input_file,
                             out_file=out_path,
                             charge_method=charge_method,
                             spin=multiplicity,
                             net_charge=net_charge,
                             atom_type=atom_type,)

        # 4. clean up
        fs.clean_temp_file_n_dir([
            input_file,
            temp_dir
        ])
        return out_path

    def _calculate_charge_with(self, charge_method: str, ncaa: NonCanonicalBase,
                               cluster_job_config: Dict) -> str:
        """generate the charge file of {charge_method} for antechamber_ncaa_to_moldesc.
        return the output file containing the charge information varied by charge_method."""
        temp_chg_file = fs.get_valid_temp_name(
            f"{eh_config['system.SCRATCH_DIR']}/temp_charge.out")
        temp_opt_file = fs.get_valid_temp_name(
            f"{eh_config['system.SCRATCH_DIR']}/temp_opt.out")
        supported_list = ["resp"]
        if charge_method in supported_list:
            if charge_method == "resp":
                # 1. optimize (H only) "PBE1PBE/def2SVP em=gd3"? Do we?
                gaussian_interface.gaussain_optimize(ncaa,
                    out_file=temp_opt_file,
                    method="PBE1PBE/def2SVP em=gd3 nosymm",
                    cluster_job_config=cluster_job_config)
                # 2. calculate charge
                gaussian_interface.gaussain_single_point(temp_opt_file,
                    out_file=temp_chg_file,
                    method="hf/6-31g* SCF=tight Pop=MK iop(6/33=2)",
                    addition_output=True,
                    cluster_job_config=cluster_job_config,)
                # 3. clean up
                fs.clean_temp_file_n_dir(temp_opt_file)
        else:
            _LOGGER.error(f"using unsupported charge_method {charge_method}.")
            raise ValueError

        return temp_chg_file

    def build_md_step(self,
                      name: str = "default",
                      # simulation
                      length: Union[float, None] = None, # ns
                      timestep: float = "default", # ns
                      minimize: bool = "default",
                      temperature: Union[float, List[Tuple[float]]] = "default",
                      thermostat: str = "default",
                      pressure_scaling: str = "default",
                      constrain: List[StructureConstraint] = "default",
                      restart: bool = "default",
                      # simulation (alternative)
                      amber_md_in_file: str = None,
                      # execution
                      core_type: str = "default",
                      cluster_job_config: Dict = "default",
                      # output
                      if_report: bool = False,
                      record_period: float = "default", # ns
                      keep_in_file: bool = False,
                      work_dir: str = "default",) -> AmberMDStep:
        """the constructor for AmberMDStep
        Args:
            name:
                the name of this md step. This will be used as the name tag in IO filenames.
            length:
                the simulation length of this step (unit: ns)
            timestep:
                the timestep of the simulation. (unit: ns)
            minimize:
                whether this md step is performing a minimization.
            temperature: 
                the temperature of the simulation. can be a list of 1d coordinates that indicate
                a changing temperature. e.g.: [(0,0), (0.5,300)] the 1st element is time in ns and
                the second is temperature at the time point.
            thermostat:
                the algorithm of the thermostat.
                options: [constant_energy, berendsen, anderson, langevin, oin, sin-respa, bussi]
            pressure_scaling:
                the scaling type for pressure
                options: [none, isotropic, anisotropic, semiisotropic]
            constrain:
                a list of StructureConstraint objects that indicates geometry constrains in the step.
            restart:
                whether restart this md from the volecity of a previous md step.
            core_type:
                the type of computing core that runs the MD. This will affect both the command
                and the cluster_job config if cluster_job_config is not None.
                options: [cpu, gpu]
            cluster_job_config:
                dictionary that assign arguments for ClusterJob.config_job
                For `res_keywords` it works as it updates the default dict in ARMerConfig.MD_GPU_RES or
                ARMerConfig.MD_CPU_RES depending on the core_type.
                NOTE that it is also used to config resources even if local run is specified.
                key list: [cluster, res_keywords]
            if_report:
                whether report result (i.e.: trajectory) of this step.
            record_period:
                if report is wanted, the simulation time period for recording a snapshot
            work_dir:
                the working dir that contains all the temp/result files.
            amber_md_in_file:
                The mdin file in Amber format. (.in) These files configures Amber execution.
                If used, settings are parsed to length,minimize,temperature,thermostat,record_period,constrain
                Note that these parsed settings can be overwritten by explicitly assign arguments.
                If `length` and `nstlim` are both specified in func argument and .in file, the value in `length`
                will be used.
            keep_in_file:
                whether the function will keep the .in file after completion.
        Return:
            AmberMDStep()
        Note:
            This function is written in a way that pro readibilty and sacrificed some ease to maintain.
            If you want to add a new argument to this function make sure you update the following sections:
                (this function)
                    0. docstring
                    1. if amber_md_in_file is not None: block
                    2. default value block
                    3. return block
                (AmberMDStep.__init__)
                    + docstring
                (AmberMDStep.md_config_dict)"""
        # pirority: direct assign > mdin file > default

        # support parse from mdin file
        if amber_md_in_file is not None:
            # parsed form {amber_md_in_file}
            # they are parsed because AmberMDStep doesn't only for generating .in files but also needs to be parsed to
            # other MolDynStep when needed. (e.g.: config OpenMM job using amber .in files)
            mdin_config = self.read_from_mdin(amber_md_in_file)
            length_mdin = mdin_config.get("length", None)
            timestep_mdin = mdin_config.get("timestep", None)
            minimize_mdin = mdin_config.get("minimize", None)
            temperature_mdin = mdin_config.get("temperature", None)
            thermostat_mdin = mdin_config.get("thermostat", None)
            pressure_scaling_mdin = mdin_config.get("pressure_scaling", None)
            record_period_mdin = mdin_config.get("record_period", None)
            constrain_mdin = mdin_config.get("constrain", None)
            restart_mdin = mdin_config.get("restart", None)

            if length is None and length_mdin is not None:
               length = length_mdin
            if minimize == "default" and minimize_mdin is not None:
               minimize = minimize_mdin
            if timestep == "default" and timestep_mdin is not None:
               timestep = timestep_mdin
            if temperature == "default" and temperature_mdin is not None:
               temperature = temperature_mdin
            if thermostat == "default" and thermostat_mdin is not None:
               thermostat = thermostat_mdin
            if pressure_scaling == "default" and pressure_scaling_mdin is not None:
               pressure_scaling = pressure_scaling_mdin
            if record_period == "default" and record_period_mdin is not None:
               record_period = record_period_mdin
            if constrain == "default" and constrain_mdin is not None:
               constrain = constrain_mdin
            if restart == "default" and restart_mdin is not None:
               restart = restart_mdin

        if length is None:
            _LOGGER.error("At least one of `length` or `nstlim` (through `amber_md_in_file`) needs to be "
                          "specified for the length of the simulation.")
            raise TypeError

        # tool: write the below code
        # print(AmberInterface._generate_default_assigning_lines_for_build_md_parameterizer_or_step(locals().items()))

        # init default values
        # make them unified in style so that amber_md_in_file block can judge if any value is explicitly assigned
        type_hint_sticker: AmberConfig
        if name == "default":
            name = self.config()["DEFAULT_MD_NAME"]
        if timestep == "default":
            timestep = self.config()["DEFAULT_MD_TIMESTEP"]
        if minimize == "default":
            minimize = self.config()["DEFAULT_MD_MINIMIZE"]
        if temperature == "default":
            temperature = self.config()["DEFAULT_MD_TEMPERATURE"]
        if thermostat == "default":
            thermostat = self.config()["DEFAULT_MD_THERMOSTAT"]
        if pressure_scaling == "default":
            pressure_scaling = self.config()["DEFAULT_MD_PRESSURE_SCALING"]
        if constrain == "default":
            constrain = self.config()["DEFAULT_MD_CONSTRAIN"]
        elif isinstance(constrain, StructureConstraint): # dispatch for non list input
            constrain = [constrain]
        if restart == "default":
            restart = self.config()["DEFAULT_MD_RESTART"]
        if core_type == "default":
            core_type = self.config()["DEFAULT_MD_CORE_TYPE"]
        if cluster_job_config == "default":
            cluster_job_config = self.config().get_default_md_cluster_job(core_type)
        else:
            # For res_keywords, it updates the default config
            cluster_job_config = copy.deepcopy(cluster_job_config)
            res_keywords_update = cluster_job_config["res_keywords"]
            default_res_keywords = self.config().get_default_md_cluster_job_res_keywords(core_type)
            cluster_job_config["res_keywords"] = default_res_keywords | res_keywords_update
        if record_period == "default":
            record_period = self.config()["DEFAULT_MD_RECORD_PERIOD_FACTOR"] * length
        if work_dir == "default":
            work_dir = self.config()["DEFAULT_MD_WORK_DIR"]
            

        return AmberMDStep(
            interface = self,
            name = name,
            length = length,
            timestep = timestep,
            minimize = minimize,
            temperature = temperature,
            thermostat = thermostat,
            pressure_scaling = pressure_scaling,
            constrain = constrain,
            restart = restart,
            core_type = core_type,
            cluster_job_config = cluster_job_config,
            if_report = if_report,
            record_period = record_period,
            keep_in_file = keep_in_file,
            work_dir = work_dir,
        )

    def convert_nc_to_mdcrd(
            self,
            nc_path: str,
            prmtop_path: str,
            out_path: str,
            autoimage: bool = True,
            remove_solvent: bool = False,
            start: int = 1,
            end: int = "last",
            step: int = 1,
            ) -> None:
        """convert a NC file to a MDCRD file using cpptraj
        Args:
            nc_path: The path to the .nc file as a str().
            prmtop_path: The path to the prmtop file as str().
            out_path: the output MDCRD file path.
            start: 1-indexed starting point. Default value is 1.
            end: 1-indexed ending point. By default uses "last"
            step: the step size of sampling points. default is 1.
            engine: The engine to convert the .nc file. Default value is "cpptraj".
        """
        contents: List[str] = [
            f"parm {prmtop_path}",
            f"trajin {nc_path} {start} {end} {step}",
        ]
        if autoimage:
            contents += [
                "autoimage"
            ]
        if remove_solvent:
            contents += [
                "strip :WAT,Cl-,Na+"
            ]
        contents += [
            f"trajout {out_path}",
            "run",
            "quit",
        ]
        contents = "\n".join(contents)
        self.run_cpptraj(contents)

    def remove_traj_solvent(
            self,
            traj_path: str,
            prmtop_path: str,
            out_path: str,
            ) -> None:
        """remove solvent in an Amber traj file using cpptraj
        Args:
            nc_path: The path to the .nc/.mdcrd file as a str().
            prmtop_path: The path to the prmtop file as str().
            out_path: the output MDCRD file path.
        """
        contents: List[str] = [
            f"parm {prmtop_path}",
            f"trajin {traj_path}",
            "strip :WAT,Cl-,Na+",
            f"trajout {out_path}",
            "run",
            "quit",
        ]
        contents = "\n".join(contents)
        self.run_cpptraj(contents)

    def load_traj(self, prmtop_path: str, traj_path: str, ref_pdb: str = None) -> StructureEnsemble:
        """load StructureEnsemble from Amber prmtop and nc/mdcrd files"""
        coord_parser_mapper = {
            ".nc" : AmberNCParser(prmtop_file=prmtop_path),
            ".mdcrd" : AmberMDCRDParser(prmtop_file=prmtop_path),
        }
        if not prmtop_io.PrmtopParser.has_add_pdb(prmtop_path):
            scratch_dir = eh_config["system.SCRATCH_DIR"]
            fs.safe_mkdir(scratch_dir)
            temp_prmtop = fs.get_valid_temp_name(f"{scratch_dir}/load_traj.prmtop")
            self.run_add_pdb(
                in_prmtop=prmtop_path,
                out_path=temp_prmtop,
                ref_pdb=ref_pdb
            )
            prmtop_path = temp_prmtop

        result = StructureEnsemble(
            topology=prmtop_path,
            top_parser=prmtop_io.PrmtopParser().get_structure,
            coordinate_list=traj_path,
            coord_parser=coord_parser_mapper[Path(traj_path).suffix].get_coordinates,
        )

        return result

    def count_num_of_frames_traj(self, prmtop_path: str, traj_path: str) -> int:
        """count the num of frames of the trajectory from
        the traj_path. Considered to be a faster approach then using
        StructureEnsemble"""
        contents: List[str] = [
            f"parm {prmtop_path}",
            f"trajin {traj_path} 1 1",
            "run",
            "quit",
        ]
        contents = "\n".join(contents)
        # prepare log
        temp_dir = eh_config.system.SCRATCH_DIR
        temp_log_path = fs.get_valid_temp_name(f"{temp_dir}/cpptraj.log")
        fs.safe_mkdir(temp_dir)
        # run
        self.run_cpptraj(contents, temp_log_path)
        # read log
        result = self.get_n_frames_cpptraj_log(temp_log_path)
        fs.clean_temp_file_n_dir([temp_log_path, temp_dir])

        return result

    # RMSD value.

    def get_rmsd(
            self,
            stru_esm: StructureEnsemble,
            stru_selection: StruSelection,
        ) -> List[float]:
        """Calculate the RMSD value of a specified StruSelection object within a StructureEnsemble instance.
        
        Args:
            stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
            stru_selection (StruSelection): A StruSelection instance representing the region for calculating RMSD value.
            
        Returns:
            rmsd_value (List[float]): A list of RMSD values for each frame in a StructureEnsemble instance comparing with the average structure.    
        """
        tmp_dir = eh_config['system.SCRATCH_DIR']
        tmp_nc_path=fs.get_valid_temp_name(os.path.join(tmp_dir, "tmp_amber_traj.nc"))
        tmp_prmtop_path=fs.get_valid_temp_name(os.path.join(tmp_dir, "tmp_amber_topology.prmtop"))
        
        self.convert_top_to_prmtop(stru_esm.topology_source_file, tmp_prmtop_path)
        self.convert_traj_to_nc(stru_esm.coordinate_list, tmp_nc_path, topology_path=tmp_prmtop_path)

        rmsd_csv_filename = fs.get_valid_temp_name(os.path.join(eh_config.system.SCRATCH_DIR, "temp_rmsd.csv"))
        amber_mask = self.get_amber_mask(stru_selection, reduce=True)
        contents: List[str] = [
            f"parm {tmp_prmtop_path}",
            f"trajin {tmp_nc_path}",
            "autoimage",
            f"rmsd {amber_mask} first mass",
            f"average crdset AVE {amber_mask}",
            "run",
            "autoimage",
            f"rmsd {amber_mask} ref AVE * out {rmsd_csv_filename} mass",
            "run",
            "quit"
        ]
        contents = "\n".join(contents)
        self.run_cpptraj(contents)

        result_df = pd.read_csv(rmsd_csv_filename, delim_whitespace=True)
        rmsd_value = result_df.iloc[:, 1]

        fs.clean_temp_file_n_dir([
            tmp_nc_path,
            tmp_prmtop_path,
            rmsd_csv_filename,
        ])
        return rmsd_value

    # -- MMPB/GBSA --
    def get_mmpbgbsa_energy(
        self,
        # Input Data
        stru_esm: StructureEnsemble,
        ligand: StruSelection,
        # Config
        work_dir: str="./binding",
        keep_in_file: bool=False,
        cluster_job_config: Union[Dict, ClusterJobConfig]= "default",
        job_check_period: int= 30, # s
        non_armer_cpu_num: int = None,
        # method specifics
        in_file: str = None,
        strip_mask: str = ":WAT,Na+,Cl-",
        solvent_model: str = "pbsa",
        igb: int = 5,
        use_sander: bool = True,
        ion_strength: float = 0.15,
        fillratio: float = 4.0,
        exdi: float = None, # amber default: 80.0
        indi: float = None, # amber default: 1.0
        ) -> List[float]:
        """Function that calculate the MMPB/GBSA energy for a given structure/structure ensemble with
        a selection pattern for the ligand.
        
        Args:
            stru, ligand, cluster_job_config, job_check_period, work_dir, keep_in_file:
                see the docstring of `enzy_htp.analysis.binding.binding_energy`
            strip_mask:
                the Amber mask for stripping unwanted parts (e.g.: solvent) in the structure ensemble
            solvent_model:
                the model for the implicit solvation. (default: pbsa)
                Supported keywords: [pbsa, gbsa]
            igb:
                the GB model to use. (default: 5 ; when solvent_model="gbsa")
                Supported value: (based on https://ambermd.org/Manuals.php)
                    0               No generalized Born term
                    1,2,5,6,7,8     Different GB models/parameters
                    3,4             Unused
                    10              numerical Poisson-Boltzmann (NOT generalized Born)
            use_sander:
                force using sander to calculate the energy
            ion_strength:
                the ion strength of the model in Molarity.
                Used as "saltcon" for gb and "istrng" for pb.
            fillratio:
                the fill ratio for PB
            in_file:
                this argument allows user to provide their own .in files. This will overwrite
                all the .in file related arguments.
        
        Returns:
            the MMPB/GBSA energies as a list"""
        # init cluster_job_config
        type_hint_sticker: AmberConfig
        default_cluster_job_config = ClusterJobConfig.from_dict(
            self.config().get_default_mmpbsa_cluster_job_config()
            )
        if cluster_job_config == "default":
            cluster_job_config = default_cluster_job_config
        elif cluster_job_config is not None:
            if not isinstance(cluster_job_config, ClusterJobConfig):
                cluster_job_config = ClusterJobConfig.from_dict(cluster_job_config)
            cluster_job_config = default_cluster_job_config | cluster_job_config
        # resolve ligand mask
        ligand_mask = self.get_amber_mask(ligand, reduce=True)

        # make MMPBSA .prmtop files
        temp_dr_prmtop = fs.get_valid_temp_name(f"{work_dir}/temp_dr.prmtop")
        temp_dl_prmtop = fs.get_valid_temp_name(f"{work_dir}/temp_dl.prmtop")
        temp_dc_prmtop = fs.get_valid_temp_name(f"{work_dir}/temp_dc.prmtop")
        temp_sc_prmtop = fs.get_valid_temp_name(f"{work_dir}/temp_sc.prmtop")
        fs.safe_mkdir(work_dir)

        self.make_mmpbgbsa_prmtop_files(
            stru_esm = stru_esm,
            ligand_mask = ligand_mask,
            strip_mask = strip_mask,
            igb = igb,
            temp_dr_prmtop = temp_dr_prmtop,
            temp_dl_prmtop = temp_dl_prmtop,
            temp_dc_prmtop = temp_dc_prmtop,
            temp_sc_prmtop = temp_sc_prmtop,
        )

        # make .nc file
        temp_nc = fs.get_valid_temp_name(f"{work_dir}/temp_mmpbsa.nc")
        self.make_mmpbgbsa_nc_file(stru_esm, temp_nc)

        # execute
        mmpbsa_result_file = fs.get_valid_temp_name(f"{work_dir}/mmpbsa_by_frames.csv")
        self.run_mmpbsa(
            dr_prmtop = temp_dr_prmtop,
            dl_prmtop = temp_dl_prmtop,
            dc_prmtop = temp_dc_prmtop,
            sc_prmtop = temp_sc_prmtop,
            traj_file = temp_nc,
            out_path = mmpbsa_result_file,
            work_dir=work_dir,
            keep_in_file = keep_in_file,
            solvent_model = solvent_model,
            cluster_job_config = cluster_job_config,
            job_check_period = job_check_period,
            non_armer_cpu_num = non_armer_cpu_num,
            in_file = in_file,
            igb = igb,
            use_sander = use_sander,
            ion_strength = ion_strength,
            fillratio = fillratio,
            exdi = exdi,
            indi = indi,
            result_in_each_frames=True,
        )

        # extract output
        result = self.parse_mmpbsa_result(mmpbsa_result_file, by_frames = True)
        result = list(result[solvent_model]["DELTA TOTAL"])

        # clean up
        fs.clean_temp_file_n_dir([
            temp_dr_prmtop,
            temp_dl_prmtop,
            temp_dc_prmtop,
            temp_sc_prmtop,
            temp_nc,
        ])

        return result

    def make_mmpbgbsa_prmtop_files(
        self,
        stru_esm: StructureEnsemble,
        ligand_mask: str,
        strip_mask: str,
        igb: int,
        temp_dr_prmtop: str,
        temp_dl_prmtop: str,
        temp_dc_prmtop: str,
        temp_sc_prmtop: str,
        ) -> None:
        """make the prmtop files needed by the MMPB/GBSA calculation and
        generate them in temp_dr_prmtop, temp_dl_prmtop, temp_dc_prmtop, temp_sc_prmtop"""
        radii = self.config()["RADII_MAP"][str(igb)]
        temp_dir = eh_config['system.SCRATCH_DIR']
        fs.safe_mkdir(temp_dir)
        temp_in_prmtop = fs.get_valid_temp_name(f"{temp_dir}/temp_in.prmtop")
        temp_in2_prmtop = fs.get_valid_temp_name(f"{temp_dir}/temp_in_2.prmtop")
        self.convert_top_to_prmtop(stru_esm.topology_source_file, temp_in_prmtop)
        self.clean_up_add_pdb_info(temp_in_prmtop, temp_in2_prmtop) # the add_pdb info in prmtop will cause bug in ante-MMPBSA
        self.run_ante_mmpbsa(
            complex_prmtop_in = temp_in2_prmtop,
            dry_complex_out = temp_dc_prmtop,
            dry_receptor_out = temp_dr_prmtop,
            dry_ligand_out = temp_dl_prmtop,
            strip_mask = strip_mask,
            radii = radii,
            ligand_mask = ligand_mask,
        )

        self.update_radii(
            prmtop_path=temp_in2_prmtop,
            out_path=temp_sc_prmtop,
            radii=radii
        )

        # clean up
        fs.clean_temp_file_n_dir([temp_in_prmtop, temp_in2_prmtop])

    def update_radii(self, prmtop_path: str, out_path: str, radii: str) -> None:
        """update the radii of the prmtop_path and generate the updated file in out_path"""

        temp_dir = eh_config['system.SCRATCH_DIR']
        fs.safe_mkdir(temp_dir)

        parmed_in_lines = [
            f"changeRadii {radii}",
            f"parmout {out_path}",
            ""
        ]
        parmed_in_str = "\n".join(parmed_in_lines)

        self.run_parmed(
            parmed_in_str,
            prmtop_path,
        )

    def make_mmpbgbsa_nc_file(
        self,
        stru_esm: StructureEnsemble,
        temp_nc: str,    
        ) -> None:
        """make the nc file for the MMPB/GBSA calculation.
        Generate the file in {temp_nc}"""

        self.convert_traj_to_nc(stru_esm.coordinate_list, temp_nc)

amber_interface = AmberInterface(None, eh_config._amber)
"""The singleton of AmberInterface() that handles all Amber related operations in EnzyHTP
Instantiated here so that other _interface subpackages can use it.
An example of this concept this AmberInterface used Gaussian for calculating the RESP charge
so it imports gaussian_interface that instantiated in the same fashion."""


class AmberRSTParser():
    """parser Amber .rst file to Structure()
    Attribute:
        prmtop_file
        parent_interface"""
    def __init__(self, prmtop_file: str, interface: BaseInterface = amber_interface):
        self.prmtop_file = prmtop_file
        self.parent_interface = interface
    
    def get_structure(self, rst_file: str) -> Structure:
        """parse a rst file to a Structure()."""


class AmberMDCRDParser():
    """parser Amber .mdcrd file
    Attribute:
        prmtop_file
        parent_interface"""
    def __init__(self, prmtop_file: str, interface: AmberInterface = amber_interface):
        self.prmtop_file = prmtop_file
        self.parent_interface = interface
    
    def get_coordinates(self, mdcrd: str, remove_solvent: bool=False) -> Generator[List[List[float]], None, None]:
        """parse a mdcrd file to a Generator of coordinates. Intermediate files are created."""
        coord = [] # coords of 1 frame
        atom_coord = [] # coord of 1 atom
        counter = 1 # collect coord every 3 count
        end_flag_1 = 0
        fake_end_flag = 0
        line_feed = os.linesep
        digit_pattern = r'[ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9]' # a number
        frame_sep_pattern = digit_pattern * 3 + line_feed

        # remove solvent
        if remove_solvent:
            temp_dir = eh_config["system.SCRATCH_DIR"]
            fs.safe_mkdir(temp_dir)
            new_mdcrd = fs.get_valid_temp_name(f"{temp_dir}/mdcrd_parser_temp_rm_sol.mdcrd")
            self.parent_interface.remove_traj_solvent(mdcrd, self.prmtop_file, new_mdcrd)
            mdcrd = new_mdcrd

        with open(mdcrd) as f:
            while True:
                # use the while True format to detect the EOF
                line=f.readline()            

                if re.match(digit_pattern, line) == None:
                    # not data line
                    if not end_flag_1 and not fake_end_flag:
                        # last line is not end line & 
                        continue
                
                if fake_end_flag:
                    # this line is an end line
                    # if next line of end is the file end
                    fake_end_flag = 0
                    if line == '':
                        break

                if end_flag_1:
                    # last line is an end line
                    if re.match(frame_sep_pattern, line) != None:
                        # this line is an end line -> last line is a fake end line
                        coord.append(holder)
                        yield coord
                        # empty for next loop
                        coord = []
                        end_flag_1 = 0
                        fake_end_flag = 1
                        continue
                    else:                        
                        # last line is a real end line
                        yield coord
                        # empty for next loop
                        coord = []
                        end_flag_1 = 0
                        if line == '':
                            #the last line
                            break
                        if line == line_feed:
                            _LOGGER.warning("unexpected empty line detected. Treat as EOF. exit reading")
                            break
                        # do not skip if normal next line

                else:
                    if re.match(frame_sep_pattern, line) != None:
                        # find line that mark the end of a frame // possible fake line
                        # store the last frame and empty the holder
                        end_flag_1 = 1 
                        lp = re.split(' +', line.strip())
                        # hold the info
                        holder = [float(lp[0]), float(lp[1]), float(lp[2])]
                        continue

                # normal data lines
                lp = re.split(' +', line.strip())
                for i in lp:
                    if counter < 3:
                        atom_coord.append(float(i))
                        counter = counter + 1
                    else:
                        atom_coord.append(float(i))
                        coord.append(atom_coord)
                        # empty for next atom
                        atom_coord = []        
                        counter = 1

    def get_structures(self, mdcrd: str) -> Generator[Structure, None, None]:
        pass

    def get_last_structure(self, mdcrd: str) -> Structure:
        """parse the last frame in a mdcrd file to a Structure(). Intermediate files are created."""
        pass


class AmberNCParser():
    """parser Amber .nc file
    Attribute:
        prmtop_file
        parent_interface"""
    def __init__(self, prmtop_file: str, interface: BaseInterface = amber_interface):
        self.prmtop_file = prmtop_file
        self.parent_interface: AmberInterface = interface

        self.mdcrd: Dict = {}
    
    def mdcrd_parser(self) -> AmberMDCRDParser:
        """get an associated mdcrd parser"""
        return AmberMDCRDParser(prmtop_file=self.prmtop_file)

    def get_coordinates(
            self,
            nc_file: str,
            autoimage: bool=True,
            remove_solvent: bool=False,
        ) -> Generator[List[List[float]], None, None]:
        """parse a nc file to a Generator of coordinates. Intermediate mdcrd file is created."""
        # 0. init temp path
        if self.mdcrd.get((autoimage, remove_solvent), None) is None:
            temp_dir = eh_config["system.SCRATCH_DIR"]
            fs.safe_mkdir(temp_dir)
            temp_mdcrd = fs.get_valid_temp_name(f"{temp_dir}/nc_parser_temp.mdcrd")

            # 1. convert to a temp mdcrd file
            self.parent_interface.convert_nc_to_mdcrd(
                nc_path=nc_file,
                prmtop_path=self.prmtop_file,
                out_path=temp_mdcrd,
                autoimage=autoimage,
                remove_solvent=remove_solvent,
            )

            # 2. store for future use
            self.mdcrd[(autoimage, remove_solvent)] = temp_mdcrd

        # 3. parse mdcrd file
        return self.mdcrd_parser().get_coordinates(self.mdcrd[(autoimage, remove_solvent)])

    def get_structures(self, nc_file: str) -> Generator[Structure, None, None]:
        pass

    def get_last_structure(self, nc_file: str) -> Structure:
        """parse the last frame in a nc file to a Structure(). Intermediate files are created."""


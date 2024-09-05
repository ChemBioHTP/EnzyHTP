"""Defines a GaussianInterface class that serves as a bridge for enzy_htp to 
utilize Gaussian software. This serves as a wrapper for all functionality provided
by the Gaussian software package. Behavior is partially controlled by the 
GaussianConfig class owned by the interface. Supported operations include:
	+ build_single_point_engine

Note that this interface is designed for Gaussian 16. The support for Gaussian 09
is uncertain. Major difference:
(new)
- opt=recalc
- MN15 functional
- geom=GIC
(default change)
- int=ACC2E=12 int=UltraFine 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-11
"""
from __future__ import annotations
import copy
import os
import re
from pathlib import Path
from typing import List, Tuple, Dict, Union
from dataclasses import dataclass
from plum import dispatch
from subprocess import CompletedProcess, SubprocessError
import goodvibes.io as gvio

from .base_interface import BaseInterface
from .handle_types import (
    QMSinglePointEngine,
    QMOptimizationEngine,
    QMResultEgg
)

from enzy_htp.core.exception import GaussianError
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.job_manager import ClusterJob
from enzy_htp.chemical import QMLevelOfTheory, MMLevelOfTheory, LevelOfTheory
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import (
    Structure,
    PDBParser,
    Residue,
    Atom,
    StructureConstraint,
    CartesianFreeze,
    StructureRegion,
    create_region_from_full_stru,
)
from enzy_htp._config.gaussian_config import GaussianConfig, default_gaussian_config
from enzy_htp import config as eh_config

@dataclass
class GaussianQMResultEgg(QMResultEgg):
    """This class define the result egg of gaussian"""
    gout_path: str
    gchk_path: str
    stru: Structure
    parent_job: ClusterJob


class GaussianChkParser():
    """parser Gaussian .chk file to wavefunctions
    Attribute:
        parent_interface"""
    def __init__(self, interface: BaseInterface):
        self.parent_interface = interface

    def get_mo_coeff(chk_file: str) -> List[List[float]]:
        """parse MO coefficients information from the chk file"""

    def get_mo_occ(chk_file: str) -> List[float]:
        """parse MO occupancy information from the chk file"""

    def get_basis_set(chk_file: str) -> str:
        """extract the basis set information from the chk file"""


class GaussianSinglePointEngine(QMSinglePointEngine):
    """The single point engine of Gaussian. QM only. Configures with
    method:
        the level of theory as QMLevelOfTheory()
    region:
        the qm region as a StructureRegion()
    keep_geom:
        whether keep/store the geometry in the result.
    cluster_job_config:
        used to make ClusterJob.
    engine:
        Gaussian16
    parent_interface:
        GaussianInterface
    name:
        the spe task name. Will be part of the file path.
    keep_in_file:
        keep gjf file or not
    work_dir:
        the working directory"""
    def __init__(self,
                interface,
                method: QMLevelOfTheory,
                region: StructureRegion,
                keep_geom: bool,
                name: str,
                cluster_job_config: Dict,
                keep_in_file: bool,
                work_dir: str,
                ):
        self._parent_interface = interface
        self._method = method
        self._region = region
        self._keep_geom = keep_geom
        self._name = name
        self._cluster_job_config = cluster_job_config
        self._keep_in_file = keep_in_file
        self._work_dir = work_dir

    # region == attribute ==
    @property
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        return "gaussian16"

    @property
    def parent_interface(self) -> GaussianInterface:
        """getter for _parent_interface"""
        return self._parent_interface

    @property
    def method(self) -> QMLevelOfTheory:
        """getter for _method"""
        return self._method

    @property
    def region(self) -> StructureRegion:
        """getter for _region"""
        return self._region

    @region.setter
    def region(self, val: StructureRegion) -> None:
        """getter for _region"""
        self._region = val

    @property
    def keep_geom(self) -> bool:
        """getter for _keep_geom"""
        return self._keep_geom

    @property
    def name(self) -> str:
        """getter for _name"""
        return self._name

    @property
    def cluster_job_config(self) -> Dict:
        """getter for _cluster_job_config"""
        return self._cluster_job_config

    @property
    def keep_in_file(self) -> bool:
        """getter for _keep_in_file"""
        return self._keep_in_file

    @property
    def work_dir(self) -> str:
        """getter for _work_dir"""
        return self._work_dir
    # endregion

    def make_job(self, stru: Structure) -> Tuple[ClusterJob, GaussianQMResultEgg]:
        """the method that makes a ClusterJob that runs the QM"""
        # 1. input
        if not isinstance(stru, Structure):
            _LOGGER.error("only allow Structure as `stru`")
            raise TypeError

        # 2. make .gjf file
        fs.safe_mkdir(self.work_dir)
        temp_gjf_file, gchk_path = self._make_gjf_file(stru)

        # 3. make cmd
        spe_cmd, gout_path = self.parent_interface.make_gaussian_cmd(temp_gjf_file)

        # 4. assemble ClusterJob
        cluster = self.cluster_job_config["cluster"]
        res_keywords = self.cluster_job_config["res_keywords"]
        env_settings = cluster.G16_ENV["CPU"]
        sub_script_path = fs.get_valid_temp_name(f"{self.work_dir}/submit_{self.name}.cmd")
        job = ClusterJob.config_job(
            commands = spe_cmd,
            cluster = cluster,
            env_settings = env_settings,
            res_keywords = res_keywords,
            sub_dir = "./", # because path are relative
            sub_script_path = sub_script_path
        )
        job.mimo = { # only used for translate clean up
            "temp_gin": [temp_gjf_file],
        }

        # 5. make result egg
        result_egg = GaussianQMResultEgg(
            gout_path = gout_path,
            gchk_path = gchk_path,
            stru=stru,
            parent_job = job,
        )

        return (job, result_egg)

    def run(self, stru: Structure) -> ElectronicStructure:
        """the method that runs the QM"""
        # 1. input
        if not isinstance(stru, Structure):
            _LOGGER.error("only allow Structure as `stru`")
            raise TypeError

        # 2. make .gjf file
        fs.safe_mkdir(self.work_dir)
        temp_gjf_file, gchk_path = self._make_gjf_file(stru)

        # 3. make cmd
        spe_cmd, gout_path = self.parent_interface.make_gaussian_cmd(temp_gjf_file)

        # 4. run cmd
        spe_cmd_exe = spe_cmd.split(" ")[0]
        spe_cmd_args = spe_cmd.split(" ")[1:]
        parent_interface = self.parent_interface
        this_spe_run = parent_interface.env_manager_.run_command(
            exe=spe_cmd_exe,
            args=spe_cmd_args,
        )
        self.check_spe_error(gout_path, this_spe_run)

        # 5. make result
        mo_parser = GaussianChkParser(
            interface=self.parent_interface)
        energy_0 = self.get_energy_0(gout_path)
        geometry = self.region.clone_to_geometry(stru)

        result = ElectronicStructure(
            energy_0 = energy_0,
            geometry = geometry,
            mo = gchk_path,
            mo_parser = mo_parser,
            source="gaussian16",
        )

        # 6. clean up
        clean_up_target = []
        if not self.keep_in_file:
            clean_up_target.append(temp_gjf_file)
        fs.clean_temp_file_n_dir(clean_up_target)

        return result

    def translate(self, result_egg: GaussianQMResultEgg) -> ElectronicStructure:
        """the method convert engine specific results to general output"""
        gout_path = result_egg.gout_path
        gchk_path = result_egg.gchk_path
        parent_job = result_egg.parent_job

        # error check
        self.check_spe_error(gout_path, result_egg.parent_job)

        stru = result_egg.stru
        mo_parser = GaussianChkParser(
            interface=self.parent_interface)
        energy_0 = self.get_energy_0(gout_path)
        geometry = self.region.clone_to_geometry(stru)

        # clean up
        clean_up_target = [result_egg.parent_job.job_cluster_log,
                           result_egg.parent_job.sub_script_path]
        if not self.keep_in_file:
            clean_up_target.extend(parent_job.mimo["temp_gin"])
        fs.clean_temp_file_n_dir(clean_up_target)

        return ElectronicStructure(
            energy_0 = energy_0,
            geometry = geometry,
            mo = gchk_path,
            mo_parser = mo_parser,
            source="gaussian16",
        )

    def check_spe_error(self,
                        gout_file: str,
                        stdstream_source: Union[ClusterJob,
                                                CompletedProcess,
                                                SubprocessError]):
        """a check for whether an error occurs in spe is needed everytime before generating
        a ElectronicStructure.
        Possible Gaussian error info places:
        1. stdout/stderr
        2. gout file
        The ultimate goal is to summarize each error type and give suggestions."""
        if not self.parent_interface.is_gaussian_completed(gout_file):
            # collect error info
            error_info_list = []
            # 1. stdout stderr
            # error types: TODO add comment here
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
            # 2. gout file
            # TODO finish this with a real example
            # normally will be a Error termination of Gaussian

            _LOGGER.error(f"Gaussian SPE didn't finish normally.{os.linesep}{os.linesep.join(error_info_list)}")
            raise GaussianError(error_info_list)

    def _make_gjf_file(self, stru: Structure) -> Tuple[str]:
        """make a temporary gjf file.
        *path* is based on self.work_dir and self.name.
        If the file exists, will change the filename to {self.name}_{index}.in
        The index start from 1.
        *content* is based on attributes of the instance.

        Return: the path of the temp gjf file."""
        # path
        temp_gjf_file_path = fs.get_valid_temp_name(f"{self.work_dir}/{self.name}.gjf")
        temp_chk_file_path = temp_gjf_file_path.replace(".gjf",".chk")

        # content
        res_keywords = self.cluster_job_config["res_keywords"]
        num_core = int(res_keywords["node_cores"])
        mem_per_core = int(res_keywords["mem_per_core"].rstrip('GB'))
        mem = (num_core * mem_per_core) - 1 # GB
        if mem < 1:
            _LOGGER.error(f"You need at least 1 GB memory to run a Gaussian job.")
            raise ValueError
        additional_keywords = []
        if self.keep_geom:
            additional_keywords = ["nosymm"]

        self.parent_interface.write_to_gjf(
            out_path = temp_gjf_file_path,
            gchk_path = temp_chk_file_path,
            num_core = num_core,
            mem = mem,
            methods = [self.method],
            job_type = "spe",
            additional_keywords = additional_keywords,
            stru = stru,
            stru_regions = [self.region],
            title = self.name,
        )

        return temp_gjf_file_path, temp_chk_file_path

    def get_energy_0(self, gout_file: str) -> float:
        """get ground state energy from the gaussian output file"""
        gout_data = self.parent_interface.read_from_gout_spe(gout_file)
        return gout_data["energy_0"]


class GaussianInterface(BaseInterface):
    """The interface between Gaussian16 and EnzyHTP.

    Attributes:
        config_: the configuration of the interface.
        env_manager_ : The EnvironmentManager() class which ensure all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: GaussianConfig = None) -> None:
        """Simplistic constructor that optionally takes a GaussianConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_gaussian_config)

    # region == mappers ==
    METHOD_KEYWORD_MAPPER = {
        "pbe0": "pbe1pbe",
        "def2-svp" : "def2svp",
        "def2-tzvp" : "def2tzvp",
        "def2-qzvpp" : "def2qzvpp",
    }
    """map Schrodinger eq. solving method name to gaussian keyword. mostly record special ones."""

    BASIS_SET_KEYWORD_MAPPER = {}
    """map basis set name to gaussian keyword."""

    SOLVENT_KEYWORD_MAPPER = {}
    """map solvent name to gaussian solvent name."""

    SOLV_MODEL_KEYWORD_MAPPER = {
        "SMD": "SMD",
    }
    """map solvation model name to gaussian keyword."""

    JOB_TYPE_KEYWORD_MAPPER = {
        "spe": "",
        "single_point": "",
        "optimization": "opt",
        "opt": "opt",
        "frequency": "freq",
        "freq": "freq",
    }
    """map job_type to gaussian keyword."""

    TAIL_SECTIONS = [
        "basis_set_spec",
        "pcm_solv_model",
    ]
    "sections after molecule specification in their orders. based on https://gaussian.com/input/?tabid=0 section ordering"

    def get_method_keyword_from_name(self, name: str) -> str:
        """map Schrodinger eq. solving method name to gaussian keyword.
        return the name itself if not in the mapper"""
        # em
        em_part = ""
        d3_pattern = r"-[dD]3$"
        d3bj_pattern = r"-[dD]3\(?[bB][jJ]\)?$"
        name = name.strip()
        if re.search(d3_pattern, name):
            name = re.sub(d3_pattern, "", name)
            em_part = " em=gd3"
        elif re.search(d3bj_pattern, name):
            name = re.sub(d3bj_pattern, "", name)
            em_part = " em=gd3bj"

        result = self.METHOD_KEYWORD_MAPPER.get(name, None)
        if result is None:
            result = name
        result = result + em_part

        return result

    def get_basis_set_keyword_from_name(self, name: str) -> Tuple[str, List[str]]:
        """map basis set name to gaussian keyword.
        return the name itself if not in the mapper
        https://gaussian.com/basissets/
        http://sobereva.com/60
        https://www.basissetexchange.org/
        Returns:
            the route component
            the tailing Gen section"""
        # TODO figure out the design of using gen and genecp when needed
        # we can use the @filepath function
        gen_section_lines = []

        result = self.BASIS_SET_KEYWORD_MAPPER.get(name, None)
        if result is None:
            result = name
        result = " "+result

        return result, gen_section_lines

    def get_solvation_keyword_from_name(self, solvent: Union[str, None], method: Union[str, None]) -> Tuple[str, List[str]]:
        """map solvation setting to gaussian keyword
        http://gaussian.com/scrf/?tabid=7
        http://sobereva.com/327"""
        # TODO figure out the design of using scrf=read probably allow mix solvent or
        # gaussian non-support solvent
        # probably fetch the supported solvent list
        # we can use the @filepath function?
        sol_read_line = []

        if solvent is None:
            result = ""
        else:
            result_sol = self.SOLVENT_KEYWORD_MAPPER.get(solvent, None)
            if result_sol is None:
                result_sol = solvent
            if method is None:
                method = self.config()["DEFAULT_SOLVENT_MODEL"]
                _LOGGER.info(f"using default solvent model: {method}")
            result_sol_method = self.SOLV_MODEL_KEYWORD_MAPPER.get(method, None)
            if result_sol_method is None:
                result_sol_method = method
            result = f" scrf=({result_sol_method},solvent={result_sol})"

        return result, sol_read_line

    # endregion

    # region == general Gaussian app interface ==
    # -- g16 --
    def get_gaussian_executable(self) -> str:
        return self.config()["GAUSSIAN_EXE"]

    def make_gaussian_cmd(self, temp_gjf_file_path) -> Tuple[str]:
        """compile the g16 cmd from config. This is general for all gaussian
        job types.
        TODO add more options when needed."""
        executable = self.get_gaussian_executable()
        temp_gout_file_path = temp_gjf_file_path.replace(".gjf", ".out")  # overwrite existing since gjf have unique name
        cmd = (f"{executable} < {temp_gjf_file_path} > {temp_gout_file_path}")
        return cmd, temp_gout_file_path

    def write_to_gjf(
        self,
        out_path: str,
        methods: List[LevelOfTheory],  # route too
        # mol spec
        stru: Structure,
        stru_regions: List[StructureRegion],
        constraints: Union[List[StructureConstraint], None] = None,
        # link0
        num_core: Union[str, None] = None,
        mem: Union[str, None] = None,
        gchk_path: Union[str, None] = None,
        # route
        job_type: Union[str, None] = None,
        job_type_additional_keyword: Union[List[str], None] = None,
        additional_keywords: Union[List[str], None] = None,
        geom_all_check: bool = False,
        # title
        title: str = "Title",
        # TODO add more when needed.
        # be the interface between actual syntax and general concept
    ) -> str:
        """function for writing the gaussian input file .gjf
        format reference: https://gaussian.com/input/?tabid=0
        This function is used in a general manner for SPE, OPT, and QMMM.
        Args:
            out_path:
                the outputing gjf path.
            method:
                the level of theory.
            stru:
                the target molecule described by Structure().
            stru_region:
                the QM region of the target molecule described by StructureRegion().
            num_core:
                the number of cpu used
            mem:
                the total memory used
            gchk_path:
                the path of the chk file
            job_type:
                lowercase job type name
                defined here: https://gaussian.com/capabilities/?tabid=1
            job_type_additional_keyword:
                job type specific keywords that dont have a good general name
                for example: the stepsize for optimization. 
            additional_keywords:
                place holder for keywords that dont have a good general name
            title:
                the title
        Return:
            (write content to out_path)
            the gjf file path"""
        gjf_lines = []
        tail_sections = {}

        # link 0
        if num_core:
            gjf_lines.append(f"%nprocshared={num_core}")
        if mem:
            gjf_lines.append(f"%mem={mem}GB")
        if gchk_path:
            gjf_lines.append(f"%chk={gchk_path}")

        # route
        if job_type_additional_keyword is None:
            job_type_additional_keyword = []
        if additional_keywords is None:
            additional_keywords = []

        route_line, tail_sections_from_route = self._make_route_line(
            methods=methods,
            job_type=job_type,
            job_type_additional_keyword=job_type_additional_keyword,
            additional_keywords=additional_keywords,
            geom_all_check=geom_all_check,
        )
        gjf_lines.extend([
            route_line,
            "",  # need file blank line
        ])
        tail_sections = tail_sections | tail_sections_from_route

        if constraints is None:
            constraints = []

        if not geom_all_check:
            # title
            gjf_lines.extend([
                title,
                "",  # need file blank line
            ])

            # mol spec
            mol_spec_lines = self._make_mol_spec(
                stru=stru,
                stru_regions=stru_regions,
                constraints=constraints,
            )
            gjf_lines.extend(mol_spec_lines)
            gjf_lines.append("")

        # TODO tailing section
        for tail_sec in self.TAIL_SECTIONS:
            sec_line = tail_sections.get(tail_sec, None)
            if sec_line is not None:
                gjf_lines.extend(sec_line)

        fs.write_lines(out_path, gjf_lines)

        return out_path

    def _make_route_line(
        self,
        methods: List[LevelOfTheory],
        job_type: str,
        job_type_additional_keyword: List[str],
        additional_keywords: str,
        geom_all_check: bool,
    ) -> Tuple[str, Dict[str, List]]:
        """make the route line for args. Only used in write_to_gjf()
        Return:
            route_line,
            tail_sections,"""
        if len(methods) > 1:
            _LOGGER.error("dont support QMMM & specifying multiple method yet. TODO")
            raise Exception("TODO")
        else:
            # lot_kw
            method = methods[0]
            lot_kw, bs_gen_lines, sol_read_line = self.lot_to_keyword(method)
            # job_type_kw
            if job_type in self.JOB_TYPE_KEYWORD_MAPPER:
                job_type_kw = self.JOB_TYPE_KEYWORD_MAPPER[job_type]
            else:
                _LOGGER.error(f"specified job_type {job_type} is not supported. Supported list: {self.JOB_TYPE_KEYWORD_MAPPER.keys()}")
                raise ValueError
            add_jobtype_kw_section = ""
            if job_type_additional_keyword:
                add_jobtype_kw_section = ",".join(job_type_additional_keyword)
                add_jobtype_kw_section = f"=({add_jobtype_kw_section})"
            # add_kw
            add_kw = " ".join(additional_keywords)
            # result
            route_line = f"#{job_type_kw}{add_jobtype_kw_section} {lot_kw} {add_kw}"
            if geom_all_check:
                route_line = f"{route_line} geom=allcheck"

            tail_sections = {
                "basis_set_spec": bs_gen_lines,
                "pcm_solv_model": sol_read_line,
            }

            return route_line, tail_sections

    def lot_to_keyword(self, lot: LevelOfTheory) -> Tuple[str, List[str], List[str]]:
        """convert a single level of theory object to gaussian keywords"""
        bs_gen_lines = []
        sol_read_line = []
        if lot.lot_type == "mm":
            raise Exception("TODO")
        else:
            lot: QMLevelOfTheory
            method_kw = self.get_method_keyword_from_name(lot.method)
            bs_kw, bs_gen_lines = self.get_basis_set_keyword_from_name(lot.basis_set)
            sol_kw, sol_read_line = self.get_solvation_keyword_from_name(lot.solvent, lot.solv_method)
            lot_kw = f"{method_kw}{bs_kw}{sol_kw}"

        return lot_kw, bs_gen_lines, sol_read_line

    def _make_mol_spec(self, stru: Structure, stru_regions: StructureRegion, constraints: StructureConstraint) -> List[str]:
        """only cartesian freeze from constraints is relevent in this part"""
        mol_spec_lines = []
        if len(stru_regions) > 1:
            _LOGGER.error("dont support QMMM & specifying multiple region yet. TODO")
            raise Exception("TODO")
        else:
            stru_region = stru_regions[0]
            if stru_region is None:
                stru_region = create_region_from_full_stru(stru)
            # chrg spin
            charge = stru_region.get_net_charge()
            spin = stru_region.get_spin()
            mol_spec_lines.append(f"{charge} {spin}")
            # geom
            # deal with constraint
            cart_freeze = []
            for cons in constraints:
                if cons.is_cartesian_freeze:
                    cart_freeze.append(cons)
            atoms = stru_region.atoms_from_geom(stru)
            geom_lines = self.get_geom_lines(atoms, cart_freeze)
            mol_spec_lines.extend(geom_lines)

        # final new line
        mol_spec_lines.append("")

        return mol_spec_lines

    def get_geom_lines(self, atoms: List[Atom], cart_freeze: List[CartesianFreeze]) -> List[str]:
        """get the geometry lines of the mol spec section
        from a list of atoms. add -1 if it is freeezed"""
        geom_lines = []
        if cart_freeze:
            raise Exception("TODO")
            cons = merge_cartesian_freeze(cart_freeze)
            for atom in atoms:
                if cons.has_atom(atom, geom=atom.root):
                    atom_line = f"{atom.element:<5} -1 {atom.coord[0]:>15.8f} {atom.coord[1]:>15.8f} {atom.coord[2]:>15.8f}"
                else:
                    atom_line = f"{atom.element:<5} {atom.coord[0]:>15.8f} {atom.coord[1]:>15.8f} {atom.coord[2]:>15.8f}"
                geom_lines.append(atom_line)

        for atom in atoms:  # TODO also support QMMM format
            atom_line = f"{atom.element:<5} {atom.coord[0]:>15.8f} {atom.coord[1]:>15.8f} {atom.coord[2]:>15.8f}"
            geom_lines.append(atom_line)

        return geom_lines

    def is_gaussian_completed(self, gout_file: str) -> bool:
        """determine whether a Gaussian Job is terminated normally.
        1. check for 'Normal termination'"""
        if os.path.exists(gout_file):
            with open(gout_file) as f:
                return "Normal termination of Gaussian" in f.read()
        else:
            return False

    def read_from_gout_spe(self, gout: str) -> Dict:
        """read a data dictionary from a gaussian output file for a
        single point based calculation. This function is placed here
        since it will be used by other read_from_gout_xxx functions.
        
        NOTE: use GoodVibes 3.2 for now.
        
        Returns;
        {   "energy_0" : ...,
            "charge" : ..., 
            "multiplicity" : ..., }"""

        (sp_energy,
         program, 
         version_program, 
         solvation_model, 
         file, 
         charge, 
         empirical_dispersion, 
         multiplicity) = gvio.parse_data(gout)
        
        result = {
            "energy_0" : sp_energy,
            "charge" : charge, 
            "multiplicity" : multiplicity, 
        }

        return result

    def read_from_gin(self, gin: str) -> Dict:
        """read a data dictionary from a gaussian input file"""
        raise Exception("TODO")

    # -- formchk --
    def run_formchk(self, chk_file: str, out_path: str) -> str:
        """interface for running formchk
        return the out_path"""
        # san check
        fs.check_file_exists(chk_file, exit_script=False)
        if fs.get_file_ext(chk_file) != ".chk":
            _LOGGER.error(f"formchk expect .chk file. (got: {chk_file})")
            raise ValueError

        self.env_manager_.run_command(self.config_.FORMCHK_EXE, [chk_file, out_path])
        return out_path

    # -- cubegen --
    def run_cubegen(self):
        """interface for running cubegen"""
        raise Exception("TODO")

    # endregion

    # region == engines ==
    def build_single_point_engine(
        self,
        # calculation config
        method: QMLevelOfTheory = "default",
        region: Union[StructureRegion, None] = None,
        keep_geom: bool = "default",
        # calculation config (alternative)
        gjf_file: Union[str, None] = None,
        # execution config
        name: str = "default",
        cluster_job_config: Dict = "default",
        keep_in_file: bool = False,
        work_dir: str = "default",
    ) -> GaussianSinglePointEngine:
        """constructor for GaussianSinglePointEngine(). Config everything for a single
        point calculation besides the actual geometry.
        Args:
            method:
                the level of theory as QMLevelOfTheory()
            region:
                the qm region as a StructureRegion()
            keep_geom:
                whether keep/store the geometry in the result.
            gjf_file:
                The input file in Gaussian format. (.gjf) These files configures Gaussian execution.
                If used, settings are parsed to
                    method
                    region
                    keep_geom
                Note that these parsed settings can be overwritten by explicitly assign arguments.
            cluster_job_config:
                dictionary that assign arguments for ClusterJob.config_job
                For `res_keywords` it works as it updates the default dict in ARMerConfig.MD_GPU_RES or
                ARMerConfig.MD_CPU_RES depending on the core_type.
                NOTE that it is also used to config resources (num cores, core type) even if local run is specified.
                key list: [cluster, res_keywords]
            keep_in_file:
                whether the function will keep the .gjf file after completion.
            work_dir:
                the working dir that contains all the temp/result files.
        Return:
            GaussianSinglePointEngine()"""
        # TODO support parse from gjf file
        if gjf_file is not None:
            # parsed form {gjf_file}
            # they are parsed because GaussianInterface doesn't only use it for generating .gjf files
            # but also needs to be parsed to other QMEngine when needed.
            # (e.g.: config XTB job using Gaussian .gjf files)
            gjf_config = self.read_from_gin(gjf_file)
            method_gjf = gjf_config.get("method", None)
            region_gjf = gjf_config.get("region", None)
            keep_geom_gjf = gjf_config.get("keep_geom", None)

            # only use when not explicitly specified
            if method == "default" and method_gjf is not None:
                method = method_gjf
            if region is None and region_gjf is not None:
                region = region_gjf
            if keep_geom == "default" and keep_geom_gjf is not None:
                keep_geom = keep_geom_gjf

        # init default values
        type_hint_sticker: GaussianConfig
        if name == "default":
            name = self.config()["DEFAULT_SPE_NAME"]
        if method == "default":
            method = self.config()["DEFAULT_SPE_METHOD"]
        if keep_geom == "default":
            keep_geom = self.config()["DEFAULT_SPE_KEEP_GEOM"]
        if cluster_job_config == "default":
            cluster_job_config = self.config().get_default_qm_spe_cluster_job_config()
        else:
            # For res_keywords, it updates the default config
            cluster_job_config = copy.deepcopy(cluster_job_config) # because we will change it in place.
            res_keywords_update = cluster_job_config["res_keywords"]
            default_res_keywords = self.config().get_default_qm_spe_cluster_job_res_keywords()
            cluster_job_config["res_keywords"] = default_res_keywords | res_keywords_update
        if work_dir == "default":
            work_dir = self.config()["DEFAULT_SPE_WORK_DIR"]

        return GaussianSinglePointEngine(
            interface=self,
            method=method,
            region=region,
            keep_geom=keep_geom,
            name=name,
            cluster_job_config=cluster_job_config,
            keep_in_file=keep_in_file,
            work_dir=work_dir,
        )

    def build_qmmm_single_point_engine():
        raise Exception("TODO")

    # endregion

    # region == TODO ==
    @dispatch
    def gaussain_optimize(self, stru: str, *args, **kwargs) -> str:
        """dispatch for using gout file as input"""
        # san check
        supported_format = [".log", ".out"]
        if not Path(stru).exists():
            _LOGGER.error(f"file dont exist: {stru}")
            raise ValueError
        if fs.get_file_ext(stru) not in supported_format:
            _LOGGER.error(f"file type not support: {stru} (supported: {supported_format})")
            raise ValueError

        # stru = GOUTParser.get_structure(stru, frame="last")
        raise Exception("TODO")

    @dispatch
    def gaussain_optimize(
        self,
        stru: Union[Residue, Structure],
        out_file: str,
        method: str,
        cluster_job_config: Dict,
        addition_output: bool = False,
    ) -> str:
        """"""
        raise Exception("TODO")

    def PDB2QMMM(
        self,
        o_dir="",
        tag="",
        work_type="spe",
        qm="g16",
        keywords="",
        prmtop_path=None,
        prepi_path: dict = None,
        spin_list=[1, 1],
        ifchk=1,
    ):
        """
        generate QMMM input template based on [connectivity, atom order, work type, layer/freeze settings, charge settings]
        * NEED TO SET LAYER BY ATOM INDEX OR SELECT A LAYER PRESET (use Config.Gaussian.layer_preset and Config.Gaussian.layer_atoms)
        * define self.frames in the func
        --------
        qm          : QM program (default: g16 / Gaussian16)
        work_type   : QMMM calculation type (default: spe)
        o_dir       : out put directory of the input file (default: self.dir/QMMM/ )
        tag         : folder name tag for potential multiple mutations
        keywords    : additional keywords add to the work_type correlated route
        prmtop_path : provide prmtop file for determining charge and spin. use self.prmtop by default
        prepi_path  : a diction of prepin file path with each ligand name as key. (e.g.: {'4CO':'./ligand/xxx.prepin'})
        spin_list   : a list of spin for each layers. (Do not support auto judge of the spin now)
        ifchk       : if save chk and return chk paths
        <see more options in Config.Gaussian>
        ========
        Gaussian
        ========
        # route section (from config module / leave work type as an inp arg)
        # charge and spin
        # coordinate
                - atom label (from .lib)
                - atom charge
                - freeze part (general option / some presupposition / freeze MM)
                - layer (same as above)
                - xyz (1. new system 2. existing template (do we still need?) -- from pdb / mdcrd / gout) ref: ONIOM_template_tool
        # connectivity
        # missing parameters (ligand related?)
        ---------
        In Config.Gaussian
        ---------
        n_cores     : Cores for gaussian job (higher pirority)
        max_core    : Per core memory in MB for gaussian job (higher pirority)
        keywords    : a list of keywords that joined together when build
        layer_preset: default preset id copied to self.layer_preset
        layer_atoms : default layer_atoms copied to self.layer_atoms
        *om_lvl     : *(can only be edit manually before loading the module) oniom method level
        """
        # san check
        support_work_type = ["spe", "opt", "tsopt"]
        if work_type not in support_work_type:
            raise Exception("PDB2QMMM.work_type : only support: " + repr(support_work_type))
        support_qm = ["g16"]
        if qm not in support_qm:
            raise Exception("PDB2QMMM.qm: only support: " + repr(support_qm))
        # default
        if prmtop_path == None:
            prmtop_path = self.prmtop_path
        # make folder
        if o_dir == "":
            o_dir = self.dir + "/QMMM" + tag
        mkdir(o_dir)
        # file path and name
        o_name = self.name + "_QMMM"
        g_temp_path = o_dir + "/" + o_name + ".gjf"
        # get stru
        self.get_stru()
        # get layer
        self._get_oniom_layer()
        # prepin path
        if prepi_path == None:
            prepi_path = self.prepi_path

        # build template
        if qm == "g16":
            self.route = self._get_oniom_g16_route(work_type, o_name, key_words=keywords)
            title = ("ONIOM input template generated by PDB2QMMM module of XXX(software name)" + line_feed)
            chrgspin = self._get_oniom_chrgspin(prmtop_path=prmtop_path, spin_list=spin_list)
            cnt_table = self.stru.get_connectivty_table(prepi_path=prepi_path)
            coord = self._get_oniom_g16_coord(prmtop_path)  # use connectivity info from the line above.
            add_prm = (self._get_oniom_g16_add_prm())  # test for rules of missing parameters

            # combine and write
            with open(g_temp_path, "w") as of:
                of.write(self.route)
                of.write(line_feed)
                of.write(title)
                of.write(line_feed)
                of.write(chrgspin)
                of.write(line_feed)
                of.write(coord)
                of.write(line_feed)
                of.write(cnt_table)
                of.write(line_feed)
                of.write(add_prm)

        # deploy to inp files
        frames = Frame.fromMDCrd(self.mdcrd)
        self.frames = frames
        gjf_paths = []
        chk_paths = []
        if Config.debug >= 1:
            print("Writing QMMM gjfs.")
        for i, frame in enumerate(frames):
            if ifchk:
                frame_path = frame.write_to_template(g_temp_path, index=str(i), ifchk=1)
                gjf_paths.append(frame_path[0])
                chk_paths.append(frame_path[1])
            else:
                gjf_paths.append(frame.write_to_template(g_temp_path, index=str(i), ifchk=0))
        # run Gaussian job
        self.qmmm_out = PDB.Run_QM(gjf_paths)

        if ifchk:
            self.qmmm_chk = chk_paths
            return self.qmmm_out, self.qmmm_chk

        return self.qmmm_out

    def _get_oniom_layer(self):
        """
        get oniom layer base on self.layer_atoms or self.layer_preset
        save a Layer object to self.layer
        """
        #  san check (need at least a set or a preset mode)
        if self.layer_atoms == [] and self.layer_preset == 0:
            raise Exception("PDB2QMMM: need layer setting. Please use self.set_oniom_layer.")
        # layer_atoms are in higher pirority
        if self.layer_atoms != []:
            self.layer = Layer(self, self.layer_atoms)
        else:
            self.layer = Layer.preset(self, self.layer_preset)

    def _get_oniom_g16_route(self, work_type, chk_name="chk_place_holder", key_words=""):
        """
        generate gaussian 16 ONIOM route section. Base on settings in the config module.
        -------
        work_type   : ONIOM calculation type (support: spe, ...)
        chk_name    : filename of chk (QMMM.chk by default / self.name + _QMMM.chk in the default workflow use.)
        key_words   : allow additional key words
        -------
        support edit keywords directly in module Config
        """
        chk = r"%chk=" + chk_name + ".chk" + line_feed
        proc = "%nprocshared=" + str(Config.n_cores) + line_feed
        mem = "%mem=" + str(Config.n_cores * Config.max_core) + "MB" + line_feed
        if type(key_words) == str and key_words != "":
            keyword_line = ("# " + " ".join(Config.Gaussian.keywords[work_type] + [
                key_words,
            ]) + line_feed)
        if type(key_words) == list:
            keyword_line = ("# " + " ".join(Config.Gaussian.keywords[work_type] + key_words) + line_feed)

        route = chk + proc + mem + keyword_line
        return route

    def _get_oniom_chrgspin(self, prmtop_path=None, spin_list=[1, 1]):
        """
        Determing charge and spin for each ONIOM layers. Base on *prmtop file* and layer settings in the *config* module.
        """
        chrgspin = None
        # san check
        if prmtop_path == None:
            if self.prmtop_path == None:
                raise Exception("Please provide or use PDB2FF() to generate a prmtop file before PDB2QMMM")
            prmtop_path = self.prmtop_path

        # get charge list
        self.chrg_list_all = PDB.get_charge_list(prmtop_path)
        # init
        self.layer_chrgspin = []
        for j in range(len(self.layer)):
            self.layer_chrgspin.append(float(0))
        # add charge to layers
        for i, chrg in enumerate(self.chrg_list_all):
            for j, layer in enumerate(self.layer):
                if i + 1 in layer:
                    self.layer_chrgspin[j] += chrg

        # add spin
        if len(self.layer_chrgspin) != len(spin_list):
            raise Exception("spin specification need to match the layer setting. e.g.: spin_list=[h_spin, l_spin]")
        for i, spin in enumerate(spin_list):
            self.layer_chrgspin[i] = (self.layer_chrgspin[i], spin)
        # make string
        if len(self.layer) == 2:
            c1 = str(round(self.layer_chrgspin[0][0]))
            s1 = str(round(self.layer_chrgspin[0][1]))
            c2 = str(round(self.layer_chrgspin[1][0]))
            s2 = str(round(self.layer_chrgspin[1][1]))
            c3 = c2
            s3 = s2
            chrgspin = " ".join([c1, s1, c2, s2, c3, s3])
        else:
            raise Exception("Only support 2 layers writing charge and spin. Update in the future")

        if Config.debug >= 1:
            print(chrgspin)

        return chrgspin

    def _get_oniom_g16_coord(self):
        """
        generate coordinate line. Base on *structure* and layer settings in the *config* module.
        Use element name as atom type for ligand atoms since they are mostly in QM regions.
        ---------------
        for a coord line:
            - element name
                - atom name (from .lib)
                - atom charge (from self.charge_list_all)
                - freeze part (general option / some presupposition)
                - xyz (from self.stru)
                - layer (general option / some presupposition)
        """
        coord = ""

        # amber default hold the chain - ligand - metal - solvent order
        a_id = 0
        for chain in self.stru.chains:
            for res in chain:
                for atom in res:
                    a_id += 1
                    # san check
                    if atom.id != a_id:
                        raise Exception("atom id error.")
                    if atom.id in self.layer[0]:
                        coord += atom.build_oniom("h", self.chrg_list_all[atom.id - 1])
                    else:
                        # consider connection
                        cnt_info = None
                        repeat_flag = 0
                        for cnt_atom in atom.connect:
                            if cnt_atom.id in self.layer[0]:
                                if repeat_flag:
                                    raise Exception("A low layer atom is connecting 2 higher layer atoms")
                                cnt_info = [
                                    "H",
                                    cnt_atom.get_pseudo_H_type(atom),
                                    cnt_atom.id,
                                ]
                                repeat_flag = 1
                        # general low layer
                        coord += atom.build_oniom("l", self.chrg_list_all[atom.id - 1], cnt_info=cnt_info)
        for lig in self.stru.ligands:
            for atom in lig:
                a_id += 1
                if atom.id != a_id:
                    raise Exception("atom id error.")
                if atom.id in self.layer[0]:
                    coord += atom.build_oniom("h", self.chrg_list_all[atom.id - 1], if_lig=1)
                else:
                    if Config.debug >= 1:
                        print("\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found ligand atom in low layer \033[0m")
                    # consider connection
                    cnt_info = None
                    repeat_flag = 0
                    for cnt_atom in atom.connect:
                        if repeat_flag:
                            raise Exception("A low layer atom is connecting 2 higher layer atoms")
                        if cnt_atom.id in self.layer[0]:
                            if Config.debug >= 1:
                                print("\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found ligand atom" + str(atom.id) +
                                      " in seperate layers \033[0m")
                            cnt_info = [
                                "H",
                                cnt_atom.get_pseudo_H_type(atom),
                                cnt_atom.id,
                            ]
                            repeat_flag = 1
                    coord += atom.build_oniom(
                        "l",
                        self.chrg_list_all[atom.id - 1],
                        cnt_info=cnt_info,
                        if_lig=1,
                    )
        for atom in self.stru.metalatoms:
            a_id += 1
            if atom.id != a_id:
                raise Exception("atom id error.")
            if atom.id in self.layer[0]:
                coord += atom.build_oniom("h", self.chrg_list_all[atom.id - 1])
            else:
                coord += atom.build_oniom("l", self.chrg_list_all[atom.id - 1])
        for sol in self.stru.solvents:
            for atom in sol:
                a_id += 1
                if atom.id != a_id:
                    raise Exception("atom id error.")
                if atom.id in self.layer[0]:
                    coord += atom.build_oniom("h", self.chrg_list_all[atom.id - 1], if_sol=1)
                else:
                    # consider connection
                    cnt_info = None  # for future update
                    repeat_flag = 0
                    for cnt_atom in atom.connect:
                        if cnt_atom.id in self.layer[0]:
                            if Config.debug >= 1:
                                print("\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found solvent atom" + str(atom.id) +
                                      " in seperate layers \033[0m")
                            if repeat_flag:
                                raise Exception("A low layer atom is connecting 2 higher layer atoms")
                            cnt_info = [
                                "H",
                                cnt_atom.get_pseudo_H_type(atom),
                                cnt_atom.id,
                            ]
                            repeat_flag = 1
                    coord += atom.build_oniom(
                        "l",
                        self.chrg_list_all[atom.id - 1],
                        cnt_info=cnt_info,
                        if_sol=1,
                    )

        return coord

    def _get_oniom_g16_add_prm(self):
        """
        Add missing parameters for protein and custom atom types
        1. addition parameters for metal element that not exist in ff96
        2. Commonly missing line for no reason: 'HrmBnd1    N   CT   HC     35.0000     109.5000'
        3. What if ligand or artificial residue appears in the low layer TODO
        4. missing parameters brought by the pseudo boundary H
        """
        # 1
        add_prm = "HrmBnd1    N   CT   HC     35.0000     109.5000" + line_feed
        # 2
        atom_rec = []
        for atom in self.stru.metalatoms:
            if atom.parm != None:
                if atom.parm[0] not in atom_rec:
                    add_prm += "VDW   " + "   ".join(atom.parm) + line_feed
                    atom_rec.append(atom.parm[0])
        # 3
        # TODO
        # 4

        return add_prm

    # endregion

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass


gaussian_interface = GaussianInterface(None, eh_config._gaussian)
"""The singleton of GaussianInterface() that handles all Gaussian related operations in EnzyHTP
Instantiated here so that other _interface subpackages can use it.
An example of this concept this AmberInterface used Gaussian for calculating the RESP charge
so it imports gaussian_interface from here."""

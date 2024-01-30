"""Defines a MultiwfnInterface class that serves as a bridge for enzy_thp to utilize the Multiwfin software. This 
serves as a wrapper for all associated Multiwfin functionality though this behavior is also partially controlled
by the MultiwfnConfig class owned by the interface. Supported operations include:
    + get_bond_dipole (bond dipole calculations)

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-01
"""
from collections import defaultdict
import re
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Union
import numpy as np

from .base_interface import BaseInterface
from .gaussian_interface import gaussian_interface

from enzy_htp.core import file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.job_manager import ClusterJob
from enzy_htp._config.multiwfn_config import MultiwfnConfig, default_multiwfn_config
from enzy_htp import config as eh_config
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import Atom


class MultiwfnInterface(BaseInterface):
    """
    Attributes:
        env_manager_:
        compatible_env_:
    """

    def __init__(self, parent, config: MultiwfnConfig = None) -> None:
        """Simplistic constructor that optionally takes an MultiwfnConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_multiwfn_config)

    # region == general Multiwfn app interface ==
    def get_multiwfn_executable(self) -> str:
        return self.config()["EXE"]

    def run_multiwfn(
        self,
        wfn_file: str,
        instr_file: str,
        log_path: str,
        cluster_job_config: Union[None, Dict],
        job_check_period: int = None,
        additional_settings: List[str] = None,
    ) -> Union[None, ClusterJob]:
        """execute Multiwfn. Also handle the ARMer or not dispatch
        Args:
            wfn_file
                the wavefunction file
            instr_file
                the instruction file
            log_path
                the stdout/err collecting file
        cluster_job_config: 
            the config for cluster_job. If None is used, the calculation
            will be run locally.
        job_check_period:
            the time cycle for update job state change (Unit: s)
        additional_settings:
            these settings will be writing to settings.ini under the work dir during
            runtime.

        Returns:
            cluster_job_config = None
                return None
            cluster_job_config != None
                return the job for running Multiwfn

        Details:
            Multiwfn {wfn_file} < {instr_file} > output
            the result location depends on the app.
            (e.g.: for bond dipole it will be in ./LMOdip.txt)"""
        # deal with non-accepting wfn_files
        temp_path_list = []
        if fs.get_file_ext(wfn_file) == ".chk":
            temp_dir = eh_config['system.SCRATCH_DIR']
            fs.safe_mkdir(temp_dir)
            temp_fchk_file = fs.get_valid_temp_name(f"{temp_dir}/temp_for_multiwfn.fchk")
            wfn_file = gaussian_interface.run_formchk(wfn_file, temp_fchk_file)
            temp_path_list.extend([temp_fchk_file, temp_dir])

        # make settings.ini

        # region deprocated NOTE we didnt use this for fchk conversion because multiwfn have a string length limit I guess so
        # it mess up with the cmd when the path is too long.
        # multiwfn_settings = "settings.ini"
        # if Path(multiwfn_settings).exists():
        #     with open(multiwfn_settings, "r+") as f:
        #         if "formchkpath" not in f.read():
        #             f.write(f'formchkpath= "{shutil.which("formchk")}"')
        # else:
        #     with open(multiwfn_settings, "w") as f:
        #         f.write(f'formchkpath= "{shutil.which("formchk")}"')
        # endregion

        if additional_settings:
            raise Exception("TODO")

        multiwfn_exe = self.get_multiwfn_executable()
        multiwfn_args = f"{wfn_file} < {instr_file} > {log_path} 2>&1"
        multiwfn_cmd = f"{multiwfn_exe} {multiwfn_args}"

        # dispatch for running types
        if cluster_job_config is None:
            # running locally
            this_spe_run = self.env_manager_.run_command(
                exe=multiwfn_exe,
                args=multiwfn_args,
            )
            self.check_multiwfn_error(log_path, this_spe_run)
            # clean up
            fs.clean_temp_file_n_dir(temp_path_list)
            return None

        elif job_check_period is None:
            _LOGGER.error("please specify job_check_period when cluster_job_config is used")
            raise ValueError
        else:
            # make job
            cluster = cluster_job_config["cluster"]
            res_keywords = cluster_job_config["res_keywords"]
            env_settings = cluster.MULTIWFN_ENV["CPU"]
            sub_script_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/submit_multiwfn.cmd")
            fs.safe_mkdir(eh_config['system.SCRATCH_DIR'])
            job = ClusterJob.config_job(
                commands=multiwfn_cmd,
                cluster=cluster,
                env_settings=env_settings,
                res_keywords=res_keywords,
                sub_dir="./",  # because path are relative
                sub_script_path=sub_script_path)
            job.mimo = {
                "temp_path_list" : temp_path_list
            }
            return job

    def check_multiwfn_error(self, out_file, stdstream_source):
        """TODO make this error checker when encounters a need."""
        pass

    # endregion

    #region engines
    def get_bond_dipole(self,
                        ele_stru: ElectronicStructure,
                        atom_1: Atom,
                        atom_2: Atom,
                        work_dir: str,
                        keep_in_file: bool,
                        cluster_job_config: Dict = "default",
                        job_check_period: int = 60,
                        **kwargs) -> Tuple[float, Tuple[float]]:
        """get bond dipole using wfn analysis with fchk files.
        TODO we also want to support single-center in some day?

        Args:
            ele_stru:
                the ElectronicStructure of a QM region containing the
                target bond.
            atom_1:
                the 1st atom defining target bond (in the region)
            atom_2:
                the 2st atom defining target bond (in the region)
            method:
                the method keyword specifying the algorithm & software
                of the dipole calculation. (see Details for supported)
            work_dir:
                the working dir that contains all the files in the SPE process
            keep_in_file:
                whether keep the input file of the calculation
            cluster_job_config: 
                the config for cluster_job. If None is used, the calculation
                will be run locally.
            job_check_period:
                the time cycle for update job state change (Unit: s)
        Returns:
            (dipole_norm_signed, dipole_vec)
                *dipole_norm_signed*
                    is the signed norm of the dipole according to its projection
                    to be bond vector.
                *dipole_vec*
                    is the vector of the dipole

                Dipole positive Direction: negative(-) to positive(+).
                Result direction: atom_1 -> atom_2

        Details:
            LMO bond dipole:
            (Method: Multiwfn manual 3.22/4.19.4)
            2-center LMO dipole is defined by the deviation of the eletronic mass center relative to the bond center.
            Dipole positive Direction: negative(-) to positive(+).
            Result direction: a1 -> a2
            For non-single bond, the final dipole is the sum of all bond dipoles between the 2 atoms.
            TODO is this really what we want for a double bond case?
            REF: Lu, T.; Chen, F., Multiwfn: A multifunctional wavefunction analyzer. J. Comput. Chem. 2012, 33 (5), 580-592.
        """
        clean_targets = []
        # init cluster_job_config
        type_hint_sticker: MultiwfnConfig
        if cluster_job_config == "default":
            cluster_job_config = self.config().get_default_bond_dipole_cluster_job_config()
        elif cluster_job_config is not None:
            # For res_keywords, it updates the default config
            res_keywords_update = cluster_job_config["res_keywords"]
            default_res_keywords = self.config().get_default_bond_dipole_res_keywords()
            cluster_job_config["res_keywords"] = default_res_keywords | res_keywords_update

        # prepare path
        wfn_file = ele_stru.mo
        fs.safe_mkdir(work_dir)
        instr_file = fs.get_valid_temp_name(f"{work_dir}/{Path(wfn_file).with_suffix('.wfnin').name}")
        log_path = fs.get_valid_temp_name(str(Path(instr_file).with_suffix('.wfnlog')))
        result_file = fs.get_valid_temp_name(str(Path(instr_file).with_suffix('.wfndip')))
        # make instr file
        settings = [
            "19",
            "-8",
            "1",
            "y",
            "q",
            "",  # final new line
        ]
        fs.write_lines(instr_file, settings)

        # run Multiwfn
        job = self.run_multiwfn(
            wfn_file=wfn_file,
            instr_file=instr_file,
            log_path=log_path,
            cluster_job_config=cluster_job_config,
            job_check_period=job_check_period,
        )
        if job is not None:
            job.submit()
            job.wait_to_end(period=job_check_period)  # wont return and do array since this is pretty fast
            self.check_multiwfn_error(log_path, job)
            clean_targets.append(job.job_cluster_log)
            clean_targets.append(job.sub_script_path)
            clean_targets.extend(job.mimo["temp_path_list"])

        # collect the ./LMOdip.txt file and clean up
        fs.safe_mv("LMOdip.txt", result_file)
        fs.clean_temp_file_n_dir(["LMOcen.txt", "new.fch", log_path] + clean_targets)
        if not keep_in_file:
            fs.clean_temp_file_n_dir([instr_file])

        # get a1->a2 vector
        bond_vec = np.array(atom_2.coord) - np.array(atom_1.coord)

        # parse the dipole result
        dipoles_data = self.parse_two_center_dp_moments(result_file)
        atom_1_id = ele_stru.geometry.get_atom_index(atom_1, indexing=1)
        atom_2_id = ele_stru.geometry.get_atom_index(atom_2, indexing=1)

        dipole_vec = np.array((0.0, 0.0, 0.0))
        target_dipoles = dipoles_data[(atom_1_id, atom_2_id)]
        if len(target_dipoles) > 1:
            _LOGGER.warning(f"found multiple LMO-based bond dipole for bond {(atom_1_id, atom_2_id)}. "
                            "The final bond dipole will be sum of all these dipole vectors. "
                            "Make sure this is what you want!")

        for this_dipole_vec in target_dipoles:
            dipole_vec += this_dipole_vec

        dipole_norm = np.linalg.norm(dipole_vec)
        if np.dot(dipole_vec, bond_vec) <= 0:
            dipole_norm *= -1

        return dipole_norm, dipole_vec

    def parse_two_center_dp_moments(self, fname: str) -> Dict[Tuple[float], List[np.array]]:
        """parse the content of LMOdip.txt file from Multiwfn LMO
        bond dipole calculation.
        
        Returns:
            a dictionary that maps bond-forming atom-index-pair to bond dipole data
            {(atom_1_id, atom_2_id) : [np.array(dp_x, dp_y, dp_z), ...], ...}
            Note that the same bond could have multiple LMO dipole as there could be
            more than the single bond"""
        # san check
        fs.check_file_exists(fname)

        lines = fs.lines_from_file(fname)

        bond_id_pattern = r'\( *([0-9]+)[A-Z][A-z]? *- *([0-9]+)[A-Z][A-z]? *\)'
        bond_data_pattern = r'X\/Y\/Z: *([0-9\.\-]+) *([0-9\.\-]+) *([0-9\.\-]+) *Norm:'

        # locate the target section
        idx = 0
        while lines[idx].strip() != "Two-center bond dipole moments (a.u.):":
            idx += 1
        idx += 1
        end: int = idx + 1
        while lines[end].find("Sum") == -1:
            end += 1

        # parse from the target section
        result = defaultdict(list)
        for ll in lines[idx:end]:
            atom_1_id, atom_2_id = re.search(bond_id_pattern, ll).groups()
            dp_x, dp_y, dp_z = re.search(bond_data_pattern, ll).groups()
            result[(int(atom_1_id), int(atom_2_id))].append(np.array((dp_x, dp_y, dp_z), dtype=np.float64))

        return result
    #endregion

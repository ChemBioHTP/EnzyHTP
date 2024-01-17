"""Defines a MultiwfnInterface class that serves as a bridge for enzy_thp to utilize the Multiwfin software. This 
serves as a wrapper for all associated Multiwfin functionality though this behavior is also partially controlled
by the MultiwfnConfig class owned by the interface. Supported operations include:
    + get_bond_dipole (bond dipole calculations)

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-01
"""
import re
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Union
import numpy as np

from .base_interface import BaseInterface

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

        Returns:
            cluster_job_config = None
                return None
            cluster_job_config != None
                return the job for running Multiwfn

        Details:
            Multiwfn {wfn_file} < {instr_file} > output
            the result location depends on the app.
            (e.g.: for bond dipole it will be in ./LMOdip.txt)"""
        multiwfn_exe = self.get_multiwfn_executable()
        multiwfn_args = f"{wfn_file} < {instr_file} > {log_path}"
        multiwfn_cmd = f"{multiwfn_exe} {multiwfn_args}"
        
        if cluster_job_config is None:
            # running locally
            this_spe_run = self.env_manager_.run_command(
                exe=multiwfn_exe,
                args=multiwfn_args,
            )
            self.check_multiwfn_error(log_path, this_spe_run)

            return None

        elif job_check_period is None:
            _LOGGER.error("please specify job_check_period when cluster_job_config is used")
            raise ValueError
        else:
            # make job
            cluster = self.cluster_job_config["cluster"]
            res_keywords = self.cluster_job_config["res_keywords"]
            env_settings = cluster.MULTIWFN_ENV["CPU"]
            sub_script_path = fs.get_valid_temp_name(
                f"{eh_config['system.SCRATCH_DIR']}/submit_multiwfn.cmd")
            job = ClusterJob.config_job(
                commands = multiwfn_cmd,
                cluster = cluster,
                env_settings = env_settings,
                res_keywords = res_keywords,
                sub_dir = "./", # because path are relative
                sub_script_path = sub_script_path
            )
            return job

    def check_multiwfn_error(self, out_file, stdstream_source):
        """TODO make this error checker when encounters a need."""
        pass
    # endregion

    #region engines
    def get_bond_dipole(
            self,
            ele_stru: ElectronicStructure,
            atom_1: Atom,
            atom_2: Atom,
            work_dir: str,
            keep_in_file: bool,
            cluster_job_config: Dict = "default",
            job_check_period: int = 60,
            **kwargs) -> List[float]:
        """get bond dipole using wfn analysis with fchk files.
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
            REF: Lu, T.; Chen, F., Multiwfn: A multifunctional wavefunction analyzer. J. Comput. Chem. 2012, 33 (5), 580-592.
        """
        
        # init cluster_job_config
        if cluster_job_config == "default":
            cluster_job_config = self.config().get_default_bond_dipole_cluster_job_config()
        elif cluster_job_config is not None:
            # For res_keywords, it updates the default config
            res_keywords_update = cluster_job_config["res_keywords"]
            default_res_keywords = self.config().get_default_qm_spe_cluster_job_res_keywords()
            cluster_job_config["res_keywords"] = default_res_keywords | res_keywords_update
        
        # prepare path
        wfn_file = ele_stru.mo
        log_path = fs.get_valid_temp_name(
            f"{work_dir}/{Path(wfn_file).with_suffix('wfnlog')}")
        result_file = fs.get_valid_temp_name(
            f"{work_dir}/{Path(wfn_file).with_suffix('wfndip')}")
        # make instr file
        instr_file = fs.get_valid_temp_name(
            f"{work_dir}/{Path(wfn_file).with_suffix('wfnin')}")
        settings = [
            "19",
            "-8",
            "1",
            "y",
            "q",
            "", # final new line
        ]
        fs.write_lines(instr_file, settings)

        # run Multiwfn
        self.run_multiwfn(
            wfn_file=wfn_file,
            instr_file=instr_file,
            log_path=log_path,
            cluster_job_config=cluster_job_config,
            job_check_period=job_check_period,
        )

        # collect the ./LMOdip.txt file and clean up
        fs.safe_mv("LMOdip.txt", result_file)
        fs.clean_temp_file_n_dir(
            ["LMOcen.txt", "new.fch", log_path]
        )
        if not keep_in_file:
            fs.clean_temp_file_n_dir(
                [instr_file]
            )

        # get a1->a2 vector
        bond_vec = np.array(atom_2.coord) - np.array(atom_1.coord)

        # parse the dipole result TODO start here
        dipoles = self.parse_two_center_dp_moments(result_file)
        for ((b1, b2), (x, y, z, norm)) in dipoles:
            if (atom_1, atom_2) != (b1, b2):
                continue
            dipole_vec = (x, y, z)
            if np.dot(np.array(dipole_vec), bond_vec) <= 0:
                norm *= -1

        return dipoles

    def parse_two_center_dp_moments(self, fname: str) -> List[Tuple[Tuple, Tuple]]:
        """parse the content of LMOdip.txt file from Multiwfn LMO
        bond dipole calculation."""
        # san check
        fs.check_file_exists(fname)

        def digits_only(raw: str) -> str:
            return re.sub(r"[a-z:/]", "", raw.lower())

        lines: List[str] = fs.lines_from_file(fname)
        idx: int = 0
        while lines[idx].find("Two-center bond dipole moments (a.u.):") == -1:
            idx += 1
        idx += 1
        end: int = idx + 1
        while lines[end].find("Sum") == -1:
            end += 1

        result = list()
        for ll in lines[idx:end]:
            tks = digits_only(ll).split()
            result.append(((int(tks[2]), int(tks[4])), tuple(map(float, tks[-4:]))))

        return result
    #endregion

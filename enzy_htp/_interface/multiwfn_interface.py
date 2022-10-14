"""Defines a MultiwfnInterface class that serves as a bridge for enzy_thp to utilize the Multiwfin software. This 
serves as a wrapper for all associated Multiwfin functionality though this behavior is also partially controlled
by the MultiwfnConfig class owned by the interface. Supported operations include:
    + bond dipole calculations
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-01
"""
import re
from pathlib import Path
from typing import List, Dict, Tuple

import numpy as np

from enzy_htp.core import env_manager as em
from enzy_htp.core import file_system as fs

# from .multiwfn_config import MultiwfnConfig, default_multiwfn_config

# TODO(CJ): add .config() getter


class MultiwfnInterface:
    """
    Attributes:
        env_manager_:
        compatible_env_:
    """

    def __init__(self, config=None):
        """ """
        self.config_ = config
        if not self.config_:
            self.config_ = default_multiwfn_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            exectubles=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()

    def parse_two_center_dp_moments(self, fname) -> List[Tuple[Tuple, Tuple]]:

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
            result.append(
                ((int(tks[2]), int(tks[4])), tuple(map(float, tks[-4:]))))

        return result

    def parse_gaussian_coords(self, fname: str):
        lines: List[str] = fs.lines_from_file(fname)
        idx: int = 0
        while lines[idx].find("Input orientation:") == -1:
            idx += 1
        idx += 5
        end: int = idx + 1
        while lines[end].find("-------------------------") == -1:
            end += 1
        atom_lines: List[str] = lines[idx:end]
        raw_coords = np.array(
            list(map(lambda ll: ll.split()[-3:], lines[idx:end])))
        return np.asfarray(raw_coords, dtype=float)

    def get_bond_dipole(self, fchks: List[str], a1, a2) -> List[float]:
        """
        get bond dipole using wfn analysis with fchk files.
        -----------
        Args:
            qm_fch_paths: paths of fchk files
                        * requires correponding out files with only ext difference
                        * (if want to compare resulting coord to original mdcrd/gjf stru)
                            requires nosymm in gaussian input that generate the fch file.
            a1          : QM I/O id of atom 1 of the target bond
            a2          : QM I/O id of atom 2 of the target bond
            prog        : program for wfn analysis (default: multiwfn)
                          **Multiwfn workflow**
                            1. Multiwfn xxx.fchk < parameter_file > output
                            the result will be in ./LMOdip.txt
                            2. extract value and project to the bond accordingly
        Returns:
            Dipoles     : A list of dipole data in a form of [(dipole_norm_signed, dipole_vec), ...]
                          *dipole_norm_signed* is the signed norm of the dipole according to its projection
                                               to be bond vector.
                          *dipole_vec* is the vector of the dipole
        -----------
        LMO bond dipole:
        (Method: Multiwfn manual 3.22/4.19.4)
        2-center LMO dipole is defined by the deviation of the eletronic mass center relative to the bond center.
        Dipole positive Direction: negative(-) to positive(+).
        Result direction: a1 -> a2
        REF: Lu, T.; Chen, F., Multiwfn: A multifunctional wavefunction analyzer. J. Comput. Chem. 2012, 33 (5), 580-592.
        """
        dipoles = []

        # self.init_Multiwfn()
        infile = f"{Path(fchks[0]).parent}/dipole_settings.in"
        settings: List[str] = (
            "19 -8 1 y q".split()
        )  # TODO(CJ): paramterize this so that it is controlled by MultiwfnConfig
        fs.write_lines(infile, settings)

        for fchk in fchks:
            # get a1->a2 vector from .out (update to using fchk TODO)
            coords = self.parse_gaussian_coords(
                str(Path(fchk).with_suffix(".out")))
            bond_vec = coords[a2] - coords[a1]

            # Run Multiwfn
            outfile = str(Path(fchk).with_suffix(".dip"))
            self.env_manager_.run_command(
                self.config_.EXE,
                [fchk, "<", infile, "&&", "mv", "LMOdip.txt", outfile])
            fs.safe_rm("LMOcen.txt")
            fs.safe_rm("new.fch")
            dp_moments = self.parse_two_center_dp_moments(outfile)
            for ((b1, b2), (x, y, z, norm)) in dp_moments:
                if (a1, a2) != (b1, b2):
                    continue
                dipole_vec = (x, y, z)
                if np.dot(np.array(dipole_vec), bond_vec) <= 0:
                    norm *= -1
            dipoles.append((norm, dipole_vec))
        return dipoles


#    @classmethod
#    def init_Multiwfn(cls, n_cores=None):
#    TODO(CJ): so this is basically a function that updates
#    the settings for the multiwfn function
#        """
#        initiate Multiwfn with settings in Config
#        """
#        # set nthreads
#        if n_cores == None:
#            n_cores = str(Config.n_cores)
#        if Config.debug >= 1:
#            print(
#                "Running: "
#                + "sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  "
#                + n_cores
#                + "/' "
#                + Config.Multiwfn.DIR
#                + "/settings.ini"
#            )
#        run(
#            "sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  "
#            + n_cores
#            + "/' "
#            + Config.Multiwfn.DIR
#            + "/settings.ini",
#            check=True,
#            text=True,
#            shell=True,
#            capture_output=True,
#        )

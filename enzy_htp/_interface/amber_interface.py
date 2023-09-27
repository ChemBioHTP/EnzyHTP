"""Defines an AmberInterface class that serves as a bridge for enzy_htp to utilize AmberMD software. Uses the AmberConfig class
found in enzy_htp/molecular_mechanics/amber_config.py. Supported operations include mutation with tLEaP, MolDynStep for module 
MD steps that can be minimization, heating, constant pressure production, or constant pressure equilibration

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-02
"""
import re
import shutil
from pathlib import Path
from subprocess import CalledProcessError
from typing import List, Tuple, Union, Dict, Any

from .base_interface import BaseInterface

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em
from enzy_htp.core.exception import UnsupportedMethod, tLEaPError
from enzy_htp._config.amber_config import AmberConfig, default_amber_config
from enzy_htp.structure.structure_io import pdb_io
import enzy_htp.structure as struct
import enzy_htp.preparation as prep
from enzy_htp import config as eh_config

class AmberInterface(BaseInterface): # TODO(qz) EOD
    """Class that provides a direct inteface for enzy_htp to utilize AmberMD software. Supported operations
    minimization, heating constant pressure production, constant pressure equilibration, trajectory file
    conversion and mutation. Users should use this class as the only way to interact with any functionality
    in Amber or associated tools like tleap.

    Attributes:
        config_	: The AmberConfig() class which provides settings for both running Amber and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensure all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: AmberConfig = None) -> None:
        """Simplistic constructor that optionally takes an AmberConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_amber_config)

    # == interface general == TODO: go to a class
    def config(self) -> AmberConfig:
        """Getter for the AmberConfig() instance belonging to the class."""
        return self.config_

    def display_config(self) -> None:
        """Prints all settings for the object's AmberConfig() inteface to stdout using AmberConfig.display()."""
        self.config_.display()

    def compatible_environment(self) -> bool:
        """Checks if the current environment is compatible with all possible needs for the AmberInterface.
        Returns:
                Whether the current environment is suitable for the AmberInterface().
        """
        return self.compatible_env_

    # == minimization-related? ==
    def write_minimize_input_file(self, fname: str, cycle: int) -> None:
        """Creates a minimization file to be used in an amber run. SHOULD NOT BE CALLED BY USERS DIRECTLY.
        All parameters in the &ctrl block are hardcoded except for ncyc and ntpr, which are 0.5*cycle
        and 0.2*cycle as integers, respectively.
        """
        minimize_lines: List[str] = [
            "Minimize",
            " &cntrl",
            "  imin=1,",
            "  ntx=1,",
            "  irest=0,",
            f"  maxcyc={str(cycle)},",
            f"  ncyc={int(0.5 * cycle)},",
            f"  ntpr={int(0.2 * cycle)},",
            "  ntwx=0,",
            "  cut=8.0,",
            " /",
            "",
        ]
        fs.write_lines(fname, minimize_lines)

    def minimize_structure(self, pdb: str, min_dir: str = "./", mode: str = "CPU", cycle: int = 2000) -> str:
        """Class method that minimizes the structure found in a supplied .pdb file, returning
        the path to the minimized file.
        Args:
            pdb: The .pdb file with the structure to be minimized.
            mode: Which version of Amber to use. Allowed values are "CPU" and "GPU".
            min_dir: Work directory where the paramter files and output are saved.
            cycle: Number of steps in the steepest descent of the cycle.
        Returns:
            Path to the minimized structure in a .pdb file.
        """
        # TODO(CJ): add some checks for inputs
        inpath = Path(pdb)
        if min_dir == "./":
            min_dir = f"{inpath.parent}/pdb_min/"

        fs.safe_mkdir(min_dir)
        outfile = f"{min_dir}/{inpath.stem}_min.pdb"
        min_in = f"{min_dir}/min.in"
        min_out = f"{min_dir}/min.out"
        min_rst = f"{min_dir}/min.ncrst"
        (prmtop, inpcrd) = self.build_param_files(str(inpath), min_dir)
        self.write_minimize_input_file(min_in, cycle)
        engine = self.config_.get_engine(mode)

        # yapf: disable
        self.env_manager_.run_command(
            engine, ["-O", "-i", min_in, "-o", min_out, "-p", prmtop, "-c", inpcrd, "-r", min_rst],

        )

        self.env_manager_.run_command(
            "ambpdb",[f"-p", prmtop, "-c", min_rst, ">", outfile]
        )
        # yapf: enable

        return outfile

    def _setup_solvation(self, igb: Union[None, str]) -> List[str]:
        """Method that creates solvation box settings for system paramterization.
        Args:
            igb: The Born solvation setting.
        Returns:
            A list() of the lines to be added to the leap.in. 
        """
        result: List[str] = list()
        if igb:
            igb_str: str = str(igb)
            radii: str = self.config_.RADII_MAP.get(igb_val)
            if radii:
                result.append(f"set default PDBRadii {radii}")
            else:
                _LOGGGER.warning(f"The IGB value {igb} is NOT suppported. Continuing...")

        result.append("center a")
        result.extend(["addions a Na+ 0", "addions a Cl- 0"])

        if self.config_.BOX_TYPE == "oct":
            result.append(f"solvateOct a TIP3PBOX {self.config_.BOX_SIZE}")
        elif self.config_.BOX_TYPE == "box":
            result.append(f"solvatebox a TIP3PBOX {self.config_.BOX_SIZE}")
        else:
            _LOGGER.error(f"The supplied box type {self.config_.BOX_SIZE} is NOT supported. Exiting...")
            exit(1)
        return result

    def build_param_files(self, in_pdb: str, build_dir: str, pH: float = 7.0, igb: str = None) -> Tuple[str, str]:
        """Creates the .prmtop and .inpcrd files for the supplied .pdb file. Handles
        processing of the Ligand() and MetalCenter() objects in the structure. Note that the supplied
        pdb file is assumed to be protonated.
        Args:
            in_pdb: The .pdb file to build parameter files for.
            buld_dir: The directory to build the parameter files in.
            pH: A float() representing the pH the parameter files should be made at. Default value is 7.0.
            igb: a str() representing Born sovlation radius settting. 
        Returns:
            A Tuple[str,str] with the containing (.prmtop path, .inpcrd path).
        """
        ligand_dir: str = f"{build_dir}/ligands/"
        metalcenter_dir: str = f"{build_dir}/metalcenter/"

        fs.safe_mkdir(ligand_dir)
        fs.safe_mkdir(metalcenter_dir)

        structure: struct.Structure = struct.PDBParser().get_structure(in_pdb)
        ligand_paths: List[str] = structure.build_ligands(ligand_dir, True)
        for lig in ligand_paths:
            _ = prep.protonate._protonate_ligand_PYBEL(lig, pH, lig)
        ligand_charges: List[int] = list(map(lambda pp: prep.protonate._ob_pdb_charge(pp), ligand_paths))
        ligand_params: List[Tuple[str, str]] = self.build_ligand_param_files(ligand_paths, ligand_charges)
        leap_path: str = f"{build_dir}/leap.in"
        leap_log: str = f"{build_dir}/leap.out"

        leap_contents: List[str] = [
            "source leaprc.protein.ff14SB",
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
        ]
        for (prepin, frcmod) in ligand_params:
            leap_contents.extend([f"loadAmberPrep {prepin}", f"loadAmberParams {frcmod}"])
        leap_contents.append(f"a = loadpdb {in_pdb}")

        leap_contents.extend(self._setup_solvation(igb))
        pdb_path: Path = Path(in_pdb)
        prmtop: str = f"{build_dir}/{pdb_path.stem}.prmtop"
        inpcrd: str = f"{build_dir}/{pdb_path.stem}.inpcrd"
        pdb_ff: str = f"{build_dir}/{pdb_path.stem}_ff.pdb"
        leap_contents.extend([f"saveamberparm a {prmtop} {inpcrd}", f"savepdb a {pdb_ff}", "quit"])
        fs.write_lines(leap_path, leap_contents)

        # yapf: disable
        self.env_manager_.run_command(
            "tleap",
            ["-s", "-f", leap_path, ">", leap_log]
        )
        # yapf: enable

        good_output: bool = True

        for pp in map(Path, [prmtop, inpcrd, pdb_ff]):
            msg = str()
            if not pp.exists():
                good_output = False
                msg = f"The file {pp} does not exist."
            else:
                if pp.stat().st_size == 0:
                    msg = f"The file {pp} exists but is empty."
                    good_output = False
            if msg:
                _LOGGER.warn(msg)

        if not good_output:
            _LOGGER.error(f"The output generated by AmberInterface.build_param_files() is incomplete. Exiting...")
            exit(1)

        return (prmtop, inpcrd)

    def build_ligand_param_files(
        self,
        paths: List[str],
        charges: List[int],
    ) -> List[Tuple[str, str]]:
        # TODO(CJ): add the method flag?
        """Creates .prepin and .frcmod files for all the supplied .pdb files. Saves files to
        same directory as the supplied .pdb files. Removes intermediate files. Should not
        be called directly by the user. Instead use AmberInterface.build_param_files()
        Args:
            paths: A list() of ligand .pdb files to prepare. MUST BE SAME LENGHT AS charges.
            charges: A list() of charges accompanying the paths. MUST BE SAME LENGTH AS paths.
        Returns:
            A list() of filename pairs with with the .prepin and .frcmod files for each supplied .pdb file.
        Raises:
            AssertionErrors: When various input sanitization checks fail.
        """
        result: List[Tuple[str, str]] = list()
        assert len(paths) == len(charges)
        for lig_pdb, net_charge in zip(paths, charges):
            lig_name: str = struct.PDBParser.get_structure(lig_pdb)
            print(lig_name)
            lig_pdb = Path(lig_pdb)
            prepin: str = str(lig_pdb.with_suffix(".prepin"))
            frcmod: str = str(lig_pdb.with_suffix(".frcmod"))
            lig_pdb: str = str(lig_pdb)
            # TODO(CJ): check if you can get the ligand name from the
            # .pdb filename alone... I think this may be possible
            # if renew
            # yapf: disable
            self.env_manager_.run_command(
                "antechamber",
                ["-i", lig_pdb, "-fi", "pdb", "-o", prepin, "-fo", "prepi", "-c", "bcc", "-s", "0", "-nc", str(net_charge), ],
            )
            # yapf: enable
            self.remove_antechamber_temp_files()
            self.env_manager_.run_command(
                # yapf: disable
                "parmchk2",
                ["-i", prepin, "-f", "prepi", "-o", frcmod])
            # yapf: enable
            # record
            result.append((prepin, frcmod))
        return result

    def md_min_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_MIN, creates a min.in file for an Amber minimization run.
        These settings are updated by accessing the owned AmberConfig object through AmberInterface.config().
        Args:
            outfile: The name of the file to save the minimization input file to..
        Returns:
            Path to the minimization file.
        """
        config = self.config_.CONF_MIN

        contents: List[str] = [
            "Minimize",
            " &cntrl",
            "  imin  = 1,",
            "  ntx   = 1,",
            "  irest = 0,",
            f"  ntc = {config['ntc']},",
            f"  ntf = {config['ntf']},",
            f"  cut   = {config['cut']},",
            f"  maxcyc = {config['maxcyc']},",
            f"  ncyc  = {int(config['ncyc_mult']*config['maxcyc'])},",
            f"  ntpr = {int(config['ntpr_mult']*config['maxcyc'])},",
            f"  ntwx = 0,",
        ]

        last_line: str = ""
        if config["ntr"] == 1:
            last_line = f"  ntr   = {config['ntr']}, restraint_wt = {config['restraint_wt']}, restraintmask = {config['restraintmask']},"

        contents.append(last_line)
        contents.extend(["/", ""])
        fs.write_lines(outfile, contents)
        return outfile

    def md_heat_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_HEAT, creates a heat.in file for an Amber heating run.
        These settings are updated by accessing the owned AmberConfig object through AmberInterface.config().
        Args:
            outfile: The name of the file to save the heating input file to.
        Returns:
            Path to the heating input file.
        """
        config = self.config_.CONF_HEAT

        contents: List[str] = [
            "Heat",
            " &cntrl",
            "  imin = 0,",
            "  ntx = 1,",
            "  irest = 0,",
            f"  ntc = {config['ntc']},",
            f"  ntf = {config['ntf']},",
            f"  cut = {config['cut']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  tempi = {config['tempi']},",
            f"  temp0 = {config['temp0']},",
            f"  ntpr =  {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx =  {int(config['ntwx']*config['nstlim'])},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 1,",
            "  ntp = 0,",
            f"  iwrap = {config['iwrap']},",
            "  nmropt = 1,",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.append(f"  ntr = {config['ntr']},")
            contents.append(f"  restraint_wt = {config['restraint_wt']},")
            contents.append(f"  restraintmask = {config['restraintmask']},")
        contents[-1] = contents[-1][0:-1]
        contents.extend([
            " /",
            " &wt",
            "  type = 'TEMP0',",
            "  istep1 = 0,",
            f"  istep2 = {int(config['A_istep2']*config['nstlim'])},",
            f"  value1 = {config['tempi']},",
            f"  value2 = {config['temp0']},",
            " /",
            f" &wt",
            "  type = 'TEMP0',",
            f"  istep1 = {int(config['A_istep2']*config['nstlim'])+1},",
            f"  istep2 = {config['nstlim']},",
            f"  value1 = {config['temp0']},",
            f"  value2 = {config['temp0']},",
            " /",
            " &wt",
            "  type = 'END',",
            " /",
            "",
        ])

        fs.write_lines(outfile, contents)
        return outfile

    def md_equi_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_EQUI, creates an input file for an Amber constant
        equilibration run. These settings are updated by accessing the owned AmberConfig object through
        AmberInterface.config().
        Args:
            outfile: The name of the file to save the constant pressure equilibration input.
        Returns:
            Path to the constant pressure equilibration input file.
        """
        config = self.config_.CONF_EQUI
        contents: List[str] = [
            "Equilibration:constant pressure",
            " &cntrl",
            "  imin = 0,",
            f"  ntx = {config['ntx']},",
            f"  irest = {config['irest']},",
            f"  ntf = {config['ntf']},",
            f"  ntc = {config['ntc']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  cut = {config['cut']},",
            f"  temp0 = {config['temp0']},",
            f"  ntpr = {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx = {config['ntwx']},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 2,",
            "  ntp = 1,",
            f"  iwrap = {config['iwrap']},",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.extend([
                f"  ntr   = {config['ntr']},",
                f"  restraint_wt = {config['restraint_wt']},",
                f"  restraintmask = {config['restraintmask']},",
            ])

        contents.extend(["/", ""])

        fs.write_lines(outfile, contents)
        return outfile

    def md_prod_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_PROD, creates an input file for an Amber production
        md run. These settings are updated by accessing the owned AmberConfig object through
        AmberInterface.config().
        Args:
            outfile: The name of the file to save the production md input.
        Returns:
            Path to the production md input file.
        """
        config = self.config_.CONF_PROD
        contents: List[str] = [
            "Production: constant pressure",
            " &cntrl",
            "  imin = 0,",
            f"  ntx = {config['ntx']},",
            f"  irest = {config['irest']},",
            f"  ntf = {config['ntf']},",
            f"  ntc = {config['ntc']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  cut = {config['cut']},",
            f"  ntpr = {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx = {config['ntwx']},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 2,",
            "  ntp = 1,",
            f"  iwrap = {config['iwrap']},",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.extend([
                f"  ntr = {config['ntr']},",
                f"  restraint_wt = {config['restraint_wt']},",
                f"  restraintmask = {config['restraintmask']}",
            ])

        contents.extend(["/", ""])
        fs.write_lines(outfile, contents)
        return outfile

    def md_run(self, prmtop: str, inpcrd: str, work_dir: str, mode: str = "CPU") -> str:
        """Runs a full MD simulation using the suplied prmtop and inpcrd files. Simulation is composed
        of four steps:
                1. minimization
                2. heating
                3. equilibration
                4. production
        The parameter files for these steps are created by AmberInterface.md_min_file(), AmberInterface.md_heat_file(),
        AmberInterface.md_equi_file() and AmberInterface.md_prod_file(), respectively. All work is done in the supplied
        work_dir location.
        Args:
            prmtop: str() with a path to a .prmtop file generated by AmberInterface.build_param_files().
            inpcrd: str() with a path to a .inpcrd file generated by AmberInterface.build_param_files().
            work_dir: Directory where the temporary files should be stored.
            mode: The mode to use for calculations. Only "CPU" and "GPU" are allowed.
        Returns:
            Path to the production .nc file.
        """
        fs.safe_mkdir(work_dir)

        min_in: str = self.md_min_file(f"{work_dir}/min.in")
        heat_in: str = self.md_heat_file(f"{work_dir}/heat.in")
        equi_in: str = self.md_equi_file(f"{work_dir}/equi.in")
        prod_in: str = self.md_prod_file(f"{work_dir}/prod.in")
        engine = self.config_.get_engine(mode)

        min_out: str = f"{work_dir}/min.out"
        min_rst: str = f"{work_dir}/min.rst"

        #yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", min_in, "-o", min_out, "-p", prmtop, "-c", inpcrd, "-r", min_rst, "-ref", inpcrd, ],
        )
        #yapf: enable

        heat_out: str = f"{work_dir}/heat.out"
        heat_rst: str = f"{work_dir}/heat.rst"

        #yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", heat_in, "-o", heat_out, "-p", prmtop, "-c", min_rst, "-ref", min_rst, "-r", heat_rst, ],
        )
        #yapf: enable

        equi_out: str = f"{work_dir}/equi.out"
        equi_rst: str = f"{work_dir}/equi.rst"
        equi_nc: str = f"{work_dir}/equi.nc"

        # yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", equi_in, "-o", equi_out, "-p", prmtop, "-c", heat_rst, "-ref", heat_rst, "-r", equi_rst, "-x", equi_nc, ],
        )
        # yapf: enable

        prod_out: str = f"{work_dir}/prod.out"
        prod_rst: str = f"{work_dir}/prod.rst"
        prod_nc: str = f"{work_dir}/prod.nc"

        # yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O",  "-i",  prod_in,  "-o",  prod_out,  "-p",  prmtop,  "-c",  equi_rst,  "-ref",  equi_rst,  "-r",  prod_rst,  "-x",  prod_nc, ],
        )
        # yapf: enable

        return prod_nc

    def run_tleap(
        self,
        tleap_in_str: str,
        # cmd based
        if_ignore_start_up: bool = True,
        additional_search_path: List[str] = None,
        tleap_out_path: str = None,
    ) -> None:
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

        NOTE: run_tleap API should not handle the index alignment since it do
        not carry information of the input pdb"""
        temp_path_list = []
        # init file paths (tleap_in_path, tleap_out_path)
        fs.safe_mkdir(eh_config["system.SCRATCH_DIR"])
        tleap_in_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/tleap.in")
        temp_path_list.extend([eh_config["system.SCRATCH_DIR"], tleap_in_path])
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

    def nc2mdcrd(
        self,
        nc_in: str,
        prmtop: str,
        point: Union[int, None] = None,
        start: int = 1,
        end: Union[int, str] = "last",
        step: int = 1,
        engine: str = "cpptraj",
    ) -> str:
        """Converts the supplied .nc file to an .mdcrd file. Scope of trajectory conversion can be
        specified though by default all available frames are selected.
        Args:
            nc_in: The path to the .nc file as a str().
            prmtop: The path to the prmtop file as  str().
            point: If not None, step is (nstlim/ntwx)/point. Default is None, but value should be int() when given.
            start: 1-indexed starting point. Default value is 1.
            end: 1-indexed ending point. By default uses "last"
            engine: The engine to convert the .nc file. Default value is "cpptraj".
        Returns:
            The name of the converted .mdcrd file.
        """
        mdcrd: str = str(Path(nc_in).with_suffix(".mdcrd"))
        config = self.config_.CONF_PROD
        if point is not None:
            step: int = int((int(config["nstlim"]) / int(config["ntwx"])) / point)

        if engine == "cpptraj":
            cpptraj_in = "./cpptraj_nc2mdcrd.in"
            cpptraj_out = "./cpptraj_nc2mdcrd.out"
            contents: List[str] = [
                f"parm {prmtop}",
                f"trajin {nc_in} {start} {end} {step}",
                f"trajout {mdcrd}",
                "run",
                "quit",
            ]
            fs.write_lines(cpptraj_in, contents)
            self.env_manager_.run_command("cpptraj", ["-i", cpptraj_in, ">", cpptraj_out])
            #fs.safe_rm(cpptraj_in)
            #fs.safe_rm(cpptraj_out)
        else:
            raise UnsupportedMethod(f"'{engine}' is not a supported method for AmberInteface.nc2mdcrd(). Allowed methods are: 'cpptraj'.")

        return mdcrd

    def get_frames(
        self,
        nc_in: str,
        prmtop: str,
        point: Union[int, None] = None,
        start: int = 1,
        end: Union[int, str] = "last",
        step: int = 1,
    ) -> str:
        """Converts the supplied .nc file to an .mdcrd file. Scope of trajectory conversion can be
        TODO(CJ): update this
        """
        outfile: str = f"{Path(prmtop).parent}/cpptraj_frames.pdb"
        config = self.config_.CONF_PROD
        if point is not None:
            step: int = int((int(config["nstlim"]) / int(config["ntwx"])) / point)

        cpptraj_in = "./cpptraj_nc2mdcrd.in"
        cpptraj_out = "./cpptraj_nc2mdcrd.out"
        contents: List[str] = [
            f"parm {prmtop}",
            f"trajin {nc_in} {start} {end} {step}",
            "strip :WAT",
            "strip :Na+,Cl-",
            f"trajout {outfile} conect",
            "run",
            "quit",
        ]
        fs.write_lines(cpptraj_in, contents)
        self.env_manager_.run_command("cpptraj", ["-i", cpptraj_in, ">", cpptraj_out])
        #fs.safe_rm(cpptraj_in)
        fs.safe_rm(cpptraj_out)
        charges = read_charge_list(prmtop)
        result = frames_from_pdb(outfile)

        if not len(result):
            return result

        num_atoms = len(result[0].atoms)
        charges = charges[0:num_atoms]
        for fidx, frame in enumerate(result):
            result[fidx].update_charges(charges)

        return result

    def parse_fmt(self, fmt: str) -> Tuple[Any, int]:
        """Private helper function that parses a %FORMAT() string from an amber .prmtop file.
        Identifies both the type and width of the format. Supplied string can contain whitespace.
        In the case that the supplied format string does not contain the pattern '%FORMAT(<>)',
        the result (None, -1) is returned.

        Args:
            fmt: the %FORMAT() line from a .prmtop as a str().

        Returns:
            A tuple() of length 2 containing the (type ctor, width as an int()).
        """

        if fmt.find('%FORMAT(') == -1 or fmt.find(')') == -1:
            return (None, -1)

        fmt = fmt.strip().replace('%FORMAT(', '').replace(')', '')

        idx = 0

        while fmt[idx].isnumeric():
            idx += 1

        dtype = None

        type_string = fmt[idx].upper()

        if type_string == 'A':
            dtype = str
        elif type_string == 'I':
            dtype = int
        elif type_string == 'E':
            dtype = float

        width = int(fmt[idx + 1:].split('.')[0])

        return (dtype, width)

    def parse_prmtop(self, prmtop: str) -> Dict:
        """Parses a .prmtop file into a dict() with the pattern of (key, value) being (keyword after %FLAG, following values). 
        Takes the specified %FORMAT() into account. Checks if the supplied file exists, exiting if it does not.

        Args:
            prmtop: str() holding the path to a .prmtop file.

        Returns:
            A dict() with (key, value) pairs of (keyword, formatted/typed values).

        """
        #TODO(CJ): unit tests

        if not Path(prmtop).exists():
            _LOGGER.error(f"The file '{prmtop}' does not exist. Exiting...")
            exit(1)

        raw: str = fs.content_from_file(prmtop)

        result = dict()

        for chunk in raw.split('%FLAG'):
            if chunk.find('%VERSION') != -1:
                #TODO(CJ): figure out what I want to do with version info
                continue

            (key, fmt, body) = chunk.split('\n', 2)
            body = re.sub('\n', '', body)
            (dtype, width) = self.parse_fmt(fmt)

            raw: List[str] = list()
            for idx in range(0, len(body), width):
                raw.append(body[idx:idx + width].strip())

            result[key.strip()] = list(map(dtype, raw))

        return result

    def remove_antechamber_temp_files(self, dname: str = './') -> None:
        """Helper method that removes temporary files generated by antechamber from a given directory.
        Removes the files "ATOMTYPE.INF", "NEWPDB.PDB", "PREP.INF" and those with the patterns "ANTECHAMBER*"
        and "sqm.*".

        Args:
            dname: str() with the path to the directory to clean out. Uses current directory as default.

        Returns:
            Nothing.

        """

        files_to_remove: List[str] = [f"{dname}/ATOMTYPE.INF", f"{dname}/PREP.INF", f"{dname}/NEWPDB.PDB"]
        files_to_remove.extend(list(map(str, Path(f"{dname}").glob("ANTECHAMBER*"))))
        files_to_remove.extend(list(map(str, Path(f"{dname}").glob("sqm*"))))

        _ = list(map(fs.safe_rm, files_to_remove))

    def add_charges(self, stru: struct.Structure, prmtop: str) -> None:
        """Method that adds RESP charges from a .prmtop file to the supplied Structure object. If the supplied prmtop
        file does not line up with the structure, an error is thrown. Performs operation in place.
        Will also throw an error if some of the resulting atoms do not get assigned charges.

        Args:
            stru: The Structure object that charges will be added to.
            prmtop: The str() path of an Amber paramter file.

        Returns:
            Nothing.
        """
        p_dict = self.parse_prmtop(prmtop)

        charges: List[float] = p_dict['CHARGE']
        anames: List[str] = p_dict['ATOM_NAME']
        res_pointers: List[int] = p_dict['RESIDUE_POINTER']
        res_names: List[int] = p_dict['RESIDUE_LABEL']

        temp = []
        for (ridx, rp), rn in zip(enumerate(res_pointers[:-1]), res_names):
            temp.extend([rn] * (res_pointers[ridx + 1] - rp))

        temp.extend((len(anames) - len(temp)) * [rn[-1]])
        res_names = temp

        for c, an, rn, aa in zip(charges, anames, res_names, stru.atoms):

            if not (an == aa.name and aa.residue.name == rn):
                _LOGGER.error(f"There is a mismatch in the atoms. Make sure the supplied prmtop file is correct. Exiting..")
                exit(1)
            aa.charge = c

        if not stru.has_charges():
            _LOGGER.error("Not all atoms in the supplied structure were assigned charges. Exiting...")
            exit(1)
            pass

    # == engines ==
    # (engines for sciencs APIs)
    def tleap_clean_up_stru(
        self,
        input_pdb: str,
        out_path: str,
        if_align_index: bool = True,
        amber_lib: str = "leaprc.protein.ff14SB",
    ) -> None:
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

"""Defines an AmberInterface class that serves as a bridge for enzy_htp to utilize AmberMD software. Uses the AmberConfig class
found in enzy_htp/molecular_mechanics/amber_config.py. Supported operations include minimization, 
heating, constant pressure production, and constant pressure equilibration.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-02
"""
from pathlib import Path
from typing import List, Tuple
import shutil
from ..core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em
from enzy_htp.core.exception import UnsupportedMethod
import enzy_htp.structure as struct
import enzy_htp.preparation as prep
from .amber_config import AmberConfig, default_amber_config


class AmberInterface:
    """Class that provides a direct inteface for enzy_htp to utilize AmberMD software.

    Atributes:
            config_	: The AmberConfig() class which provides settings.
            env_manager_ : TODO
            compatible_env_ : TODO
    """

    def __init__(self, config: AmberConfig = None) -> None:
        """Simplistic constructor that optionally takes an AmberConfig object as its only argument."""
        self.config_ = config
        if not self.config_:
            self.config_ = default_amber_config()
        self.env_manager_ = em.EnvironmentManager(
            env_vars=self.config_.required_env_vars(),
            exectubles=self.config_.required_executables(),
        )
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()

    def config(self) -> AmberConfig:
        """TODO(CJ)"""
        return self.config_

    def display_config(self) -> None:
        """TODO(CJ)"""
        self.config_.display()

    def compatible_environment(self) -> bool:
        """Checks if the current environment is compatible with all possible needs for the AmberInterface.

        Returns:
                Whether the current environment is suitable for the AmberInterface().
        """
        return self.compatible_env_

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

    def minimize_structure(
        self, pdb: str, min_dir: str = "./", mode: str = "CPU", cycle: int = 2000
    ) -> str:
        """Class method that minimizes the structure found in a supplied .pdb file, returning
        the path to the minimized file.

                Args:
                        pdb: The .pdb file with the structure to be minimized.
                        mode: Which version of Amber to use. Allowed values are "CPU" and "GPU".
                        min_dir: Work directory where the paramter files and output are saved.
                        cycle: TODO

                Returns:
                        Path to the minimized structure in a .pdb file.
        """
        # TODO(CJ): add some checks for inputs
        inpath = Path(pdb)
        min_dir = f"{inpath.parent}/pdb_min/"
        fs.safe_mkdir(min_dir)
        outfile = f"{min_dir}/{inpath.stem}_min.pdb"
        min_in = f"{min_dir}/min.in"
        min_out = f"{min_dir}/min.out"
        min_rst = f"{min_dir}/min.ncrst"
        (prmtop, inpcrd) = self.build_param_files(str(inpath), min_dir)
        self.write_minimize_input_file(min_in, cycle)
        engine = self.config_.get_engine(mode)

        self.env_manager_.run_command(
            engine,
            [
                "-O",
                "-i",
                min_in,
                "-o",
                min_out,
                "-p",
                prmtop,
                "-c",
                inpcrd,
                "-r",
                min_rst,
            ],
        )

        self.env_manager_.run_command(
            "ambpdb", [f"-p", prmtop, "-c", min_rst, ">", outfile]
        )
        # shutil.move(prmtop, min_dir )
        # shutil.move(inpcrd, min_dir )
        return outfile

    def build_param_files(self, in_pdb: str, build_dir: str) -> Tuple[str, str]:
        """Creates the .prmtop and .inpcrd files for the supplied .pdb file. Handles
        processing of the Ligand() and MetalCenter() objects in the structure.

                Args:
                        in_pdb: The .pdb file to build parameter files for.
                        buld_dir: The directory to build the parameter files in.

                Returns:
                        A Tuple[str,str] with the containing (.prmtop path, .inpcrd path).
        """
        ligand_dir: str = f"{build_dir}/ligand/"
        metalcenter_dir: str = f"{build_dir}/metalcenter/"
        fs.safe_mkdir(ligand_dir)
        fs.safe_mkdir(metalcenter_dir)
        structure: struct.Structure = struct.structure_from_pdb(in_pdb)
        ligand_paths: List[str] = structure.build_ligands(ligand_dir, True)
        # TODO(CJ): we have to protonate this first. not sure if this id documented elsewhere
        for lig in ligand_paths:
            _ = prep.protonate._protonate_ligand_PYBEL(lig, 7.0, lig)
        ligand_charges: List[int] = list(
            map(lambda pp: prep.protonate._ob_pdb_charge(pp), ligand_paths)
        )
        ligand_params: List[Tuple[str, str]] = self.build_ligand_param_files(
            ligand_paths, ligand_charges
        )
        leap_path: str = f"{build_dir}/leap.in"
        leap_log: str = f"{build_dir}/leap.out"
        # sol_path: str = f"{build_dir}/leap.in"
        leap_contents: List[str] = [
            "source leaprc.protein.ff14SB",
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
        ]
        for (prepin, frcmod) in ligand_params:
            leap_contents.extend(
                [f"loadAmberParams {frcmod}", f"loadAmberPrep {prepin}"]
            )
        leap_contents.append(f"a = loadpdb {in_pdb}")
        # TODO(CJ): Include igb
        leap_contents.append("center a")
        # TODO(CJ): include solvation stuff as a separate function
        leap_contents.extend(["addions a Na+ 0", "addions a Cl- 0"])
        if self.config_.BOX_TYPE == "oct":
            leap_contents.append(f"solvateOct a TIP3PBOX {self.config_.BOX_SIZE}")
        else:
            leap_contents.append(f"solvatebox a TIP3PBOX {self.config_.BOX_SIZE}")
        pdb_path: Path = Path(in_pdb)
        prmtop: str = f"{build_dir}/{pdb_path.stem}.prmtop"
        inpcrd: str = f"{build_dir}/{pdb_path.stem}.inpcrd"
        pdb_ff: str = f"{build_dir}/{pdb_path.stem}_ff.pdb"
        leap_contents.extend(
            [f"saveamberparm a {prmtop} {inpcrd}", f"savepdb a {pdb_ff}", "quit"]
        )
        fs.write_lines(leap_path, leap_contents)
        # TODO(CJ): Check that this actually works before returning
        self.env_manager_.run_command("tleap", ["-s", "-f", leap_path, ">", leap_log])
        return (prmtop, inpcrd)

    def build_ligand_param_files(
        self, paths: List[str], charges: List[int]
    ) -> List[Tuple[str, str]]:
        # TODO(CJ): add the method flag?
        """Creates .prepin and .frcmod files for all the supplied .pdb files. Saves files to
        same directory as the supplied .pdb files. Removes intermediate files. Should not
        be called directly by the user. Instad use TODO(CJ): name of method here

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
            lig_name: str = struct.structure_parser.get_ligand_name(lig_pdb)
            lig_pdb = Path(lig_pdb)
            prepin: str = str(lig_pdb.with_suffix(".prepin"))
            frcmod: str = str(lig_pdb.with_suffix(".frcmod"))
            lig_pdb: str = str(lig_pdb)
            # TODO(CJ): check if you can get the ligand name from the
            # .pdb filename alone... I think this may be possible
            # if renew
            # TODO(CJ): figure out how to implement the environment manager here
            self.env_manager_.run_command(
                "antechamber",
                [
                    "-i",
                    lig_pdb,
                    "-fi",
                    "pdb",
                    "-o",
                    prepin,
                    "-fo",
                    "prepi",
                    "-c",
                    "bcc",
                    "-s",
                    "0",
                    "-nc",
                    str(net_charge),
                ],
            )
            files_to_remove: List[
                str
            ] = "ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.pdb sqm.in sqm.out".split()
            files_to_remove.extend(list(map(str, Path(".").glob("ANTECHAMBER*"))))
            _ = list(map(lambda fname: fs.safe_rm(fname), files_to_remove))
            # gen frcmod
            # TODO(CJ): add some kind of check that this all actually runs correctly w/o errors
            self.env_manager_.run_command(
                "parmchk2", ["-i", prepin, "-f", "prepi", "-o", frcmod]
            )
            # record
            result.append((prepin, frcmod))
        return result

    def md_min_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_MIN, creates a min.in file for an Amber minimization run.
        TODO(CJ): add something about updating settings
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
        """TODO(CJ)"""
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
        contents.extend(
            [
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
            ]
        )

        fs.write_lines(outfile, contents)
        return outfile

    def md_equi_file(self, outfile: str) -> str:
        """TODO(CJ)"""
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
            contents.extend(
                [
                    f"  ntr   = {config['ntr']},",
                    f"  restraint_wt = {config['restraint_wt']},",
                    f"  restraintmask = {config['restraintmask']},",
                ]
            )

        contents.extend(["/", ""])

        fs.write_lines(outfile, contents)
        return outfile

    def md_prod_file(self, outfile: str) -> str:
        """TODO(CJ)"""
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
            contents.extend(
                [
                    f"  ntr = {config['ntr']},",
                    f"  restraint_wt = {config['restraint_wt']},",
                    f"  restraintmask = {config['restraintmask']}",
                ]
            )

        contents.extend(["/", ""])
        fs.write_lines(outfile, contents)
        return outfile

    # TODO(CJ): add some way to interact with the AmberConfig() config_ attribute
    def md_run(self, prmtop: str, inpcrd: str, work_dir: str, mode: str = "CPU") -> str:
        """Runs a full MD simulation using the suplied prmtop and inpcrd files. Simulation is composed
        of four steps:
                1. minimization
                2. heating
                3. equilibration
                4. production

        Args:
                prmtop:
                                inpcrd:
                                work_dir:
                                mode:

        Returns:
                Path to the

        Raises:
                TODO
        """
        fs.safe_mkdir(work_dir)
        min_in: str = self.md_min_file(f"{work_dir}/min.in")
        heat_in: str = self.md_heat_file(f"{work_dir}/heat.in")
        equi_in: str = self.md_equi_file(f"{work_dir}/equi.in")
        prod_in: str = self.md_prod_file(f"{work_dir}/prod.in")
        engine = self.config_.get_engine(mode)

        min_out: str = f"{work_dir}/min.out"
        min_rst: str = f"{work_dir}/min.rst"
        self.env_manager_.run_command(
            engine,
            [
                "-O",
                "-i",
                min_in,
                "-o",
                min_out,
                "-p",
                prmtop,
                "-c",
                inpcrd,
                "-r",
                min_rst,
                "-ref",
                inpcrd,
            ],
        )

        heat_out: str = f"{work_dir}/heat.out"
        heat_rst: str = f"{work_dir}/heat.rst"
        self.env_manager_.run_command(
            engine,
            [
                "-O",
                "-i",
                heat_in,
                "-o",
                heat_out,
                "-p",
                prmtop,
                "-c",
                min_rst,
                "-ref",
                min_rst,
                "-r",
                heat_rst,
            ],
        )

        equi_out: str = f"{work_dir}/equi.out"
        equi_rst: str = f"{work_dir}/equi.rst"
        equi_nc: str = f"{work_dir}/equi.nc"
        self.env_manager_.run_command(
            engine,
            [
                "-O",
                "-i",
                equi_in,
                "-o",
                equi_out,
                "-p",
                prmtop,
                "-c",
                heat_rst,
                "-ref",
                heat_rst,
                "-r",
                equi_rst,
                "-x",
                equi_nc,
            ],
        )

        prod_out: str = f"{work_dir}/prod.out"
        prod_rst: str = f"{work_dir}/prod.rst"
        prod_nc: str = f"{work_dir}/prod.nc"
        self.env_manager_.run_command(
            engine,
            [
                "-O",
                "-i",
                prod_in,
                "-o",
                prod_out,
                "-p",
                prmtop,
                "-c",
                equi_rst,
                "-ref",
                equi_rst,
                "-r",
                prod_rst,
                "-x",
                prod_nc,
            ],
        )

        return prod_nc

    def nc2mdcrd(
        self,
        nc_file: str,
        prmtop: str,
        point=None,
        start=1,
        end=-1,
        step=1,
        engine: str = "cpptraj",
    ) -> str:
        """Converts the supplied .nc file to an .mdcrd file.
        convert self.nc to a mdcrd file to read and operate.(self.nc[:-2]+'.mdcrd' by default)
        a easier way is to use pytraj directly.
        ---------------
        o_path: user assigned out path (self.nc[:-2]+'mdcrd' by default)
        point:  sample point. use value from self.conf_prod['nstlim'] and self.conf_prod['ntwx'] to determine step size.
        start:  start point
        end:    end point
        step:   step size
        engine: pytraj or cpptraj (some package conflict may cause pytraj not available)
        """
        mdcrd: str = str(Path(nc_file).with_suffix(".mdcrd"))
        config = self.config_.CONF_PROD
        if end == -1:
            end = "last"
        if point != None:
            all_p = int(config["nstlim"]) / int(config["ntwx"])
            step = int(all_p / point)

        if engine == "cpptraj":
            cpptraj_in = self.cache_path + "/cpptraj_nc2mdcrd.in"
            cpptraj_out = self.cache_path + "/cpptraj_nc2mdcrd.out"
            contents: List[str] = [
                f"parm {prmtop}",
                f"trajin {nc_file} {start} {end} {step}",
                f"trajout {mdcrd}",
                "run",
                "quit",
            ]
            fs.write_lines(cpptraj_in, contents)
            self.env_manager_.run_command(
                "cpptraj", ["-i", cpptraj_in, ">", cpptraj_out]
            )
        else:
            pass

        return mdcrd

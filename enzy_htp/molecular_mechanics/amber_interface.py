"""Defines an AmberInterface class that serves as a bridge for enzy_htp to utilize AmberMD software. Uses the AmberConfig class
found in enzy_htp/molecular_mechanics/amber_interface.py. Supported operations include minimization, 
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
import enzy_htp.structure as struct
import enzy_htp.preparation as prep
from .amber_config import AmberConfig, default_amber_config


class AmberInterface:
    """Class that provides a direct inteface for enzy_htp to utilize AmberMD software.

    Atributes:
            config_	: The AmberConfig class which provides settings.
    """

    def __init__(self, config: AmberConfig = None) -> None:
        """Simplistic constructor that optionally takes an AmberConfig object as its only argument."""
        self.config_ = config
        if not self.config_:
            self.config_ = default_amber_config()

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
        self, pdb: str, mode: str = "CPU", min_dir: str = "./", cycle: int = 2000
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
		#TODO(CJ): add some checks for inputs
        inpath = Path(pdb)
        min_dir = f"{inpath.parent}/pdb_min/"
        fs.safe_mkdir(min_dir)
        outfile = f"{min_dir}/{inpath.stem}_min.pdb"
        min_in = f"{min_dir}/min.in"
        min_out = f"{min_dir}/min.out"
        min_rst = f"{min_dir}/min.ncrst"
        (prmtop, inpcrd) = self.build_param_files(inpath, outdir)
        self.write_minimize_input_file(min_in, cycle)
        engine = self.config_.get_engine( mode )
        
        self.config_.env_manager().run_command(
           engine,
               ["-O", "-i", min_in, "-o", min_out, "-p", prmtop, "-c", inpcrd, "-r", min_rst,]
        )
        
        self.config_.env_manager().run_command(
           "ambpdb", [f"-p", prmtop, "-c", min_rst, ">", outfile]
        )

    
        shutil.move(prmtop, min_dir )
        shutil.move(inpcrd, min_dir )

        return outfile

    def build_param_files(self, in_pdb: str, build_dir: str) -> Tuple[str, str]:
        """Creates the .prmtop and .inpcrd files for the supplied .pdb file. Handles 
        processing of the Ligand() and MetalCenter() objects in the structure.

		Args:
			in_pdb: The .pdb file to build parameter files for.
			buld_dir: The directory to build the parameter files in.

		Returns:
			A Tuple[str,str] with the layout (.prmtop path, .inpcrd path).
		"""
        ligand_dir: str = f"{build_dir}/ligand/"
        metalcenter_dir: str = f"{build_dir}/metalcenter/"
        fs.safe_mkdir(ligand_dir)
        fs.safe_mkdir(metalcenter_dir)
        structure: struct.Structure = struct.structure_from_pdb(in_pdb)
        ligand_paths: List[str] = structure.build_ligands(ligand_dir, True)
        ligand_charges: List[int] = list(
            map(lambda pp: prep._ob_pdb_charge(pp), ligand_paths)
        )
        ligand_params: List[Tuple[str,str]] = self.build_ligand_param_files(ligand_paths,ligand_charges)
        leap_path: str = f"{build_dir}/leap.in"
        leap_log: str = f"{build_dir}/leap.out"
        #sol_path: str = f"{build_dir}/leap.in"
        leap_contents: List[str] = ['source leaprc.protein.ff14SB', 'source leaprc.gaff', 'source leaprc.water.tip3p']
        for (prepin, frcmod) in ligand_params:
            leap_contents.extend([f"loadAmberParams {frcmod}", f"loadAmberPrep {prepin}"])
        leap_contents.append(f"a = loadpdb {in_pdb}")
        #TODO(CJ): Include igb
        leap_contents.append("center a")
        #TODO(CJ): include solvation stuff
        pdb_path: Path = Path(in_pdb)
        prmtop: str = f"{build_dir}/{pdb_path.stem}.prmtop"
        inpcrd: str = f"{build_dir}/{pdb_path.stem}.inpcrd"
        pdb_ff: str = f"{build_dir}/{pdb_path.stem}_ff.pdb"
        leap_contents.extend([f"saveamberparm a {prmptop} {inpcrd}", f"savepdb a {pdb_ff}", "quit"])
        fs.write_lines(leap_path, leap_contents)
        #TODO(CJ): Check that this actually works before returning
        os.system(f"tleap -s -f {leap_path} > {leap_log}")
        return (prmtop, inpcrd)

    def PDB2FF(
        self,
        prm_out_path="",
        o_dir="",
        lig_method="AM1BCC",
        renew_lig=0,
        local_lig=1,
        ifsavepdb=0,
        igb=None,
        if_prm_only=0,
    ):
        """
        PDB2FF(self, o_dir='')
        --------------------
        prm_out_path: output path of the prmtop file
        o_dir contral where the leap.in and leap.log go: has to contain a / at the end (e.g.: ./dir/)
        renew_lig: 0 use old ligand parm files if detected.
                   1 generate new ones.
        local_lig: 0 export lig files to the workdir level in HTP jobs.
                   1 keep local
        --------------------
        chains:
        ligand: - less junk files if your workflow contains a protonation step in advance.
        metal:
        """
        # check and generate self.stru
        self.get_stru()

        # build things seperately
        if local_lig:
            lig_dir = self.dir + "/ligands/"
            met_dir = self.dir + "/metalcenters/"
        else:
            lig_dir = self.dir + "/../ligands/"
            met_dir = self.dir + "/../metalcenters/"
        mkdir(lig_dir)
        mkdir(met_dir)

        ligands_pathNchrg = self.stru.build_ligands(lig_dir, ifcharge=1, ifunique=1)
        # metalcenters_path = self.stru.build_metalcenters(met_dir)
        # parm
        ligand_parm_paths = self._ligand_parm(
            ligands_pathNchrg, method=lig_method, renew=renew_lig
        )
        # self._metal_parm(metalcenters_path)
        # combine
        if o_dir != "":
            mkdir(o_dir)
        self._combine_parm(
            ligand_parm_paths,
            prm_out_path=prm_out_path,
            o_dir=o_dir,
            ifsavepdb=ifsavepdb,
            igb=igb,
            if_prm_only=if_prm_only,
        )
        if ifsavepdb:
            self.path = self.path_name + "_ff.pdb"
            self._update_name()

        return (self.prmtop_path, self.inpcrd_path)

        def build_ligand_param_files(
            self, paths: List[str], charges: List[int]
        ) -> List[Tuple[str, str]]:
            # TODO(CJ): add the method flag?
            """TODO(CJ): Documentation"""
            result: List[Tuple[str, str]] = list()
            parm_paths = []
            self.prepi_path = {}
            # TODO(CJ): Beef up these checks
            assert len(paths) == len(charges)
            for lig_pdb, net_charge in zip(paths, charges):
                lig_pdb = Path(lig_pdb)
                prepin: Path = lig_pdb.with_suffix(".prepin")
                frcmod: Path = lig_pdb.with_suffix(".frcmod")
                # TODO(CJ): check if you can get the ligand name from the
                # .pdb filename alone... I think this may be possible
                lig_name: str = struct.structure_parser.get_ligand_name(lig_pdb)
                # if renew
                # TODO(CJ): figure out how to implement the environment manager here
                run(
                    f"{Config.Amber.AmberHome}/bin/antechamber -i {lig_pdb} -fi pdb -o {prepin} -fo prepi -c bcc -s 0 -nc {net_charge}"
                )
                files_to_remove: List[
                    str
                ] = "ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.pdb sqm.in sqm.out".split()
                files_to_remove.extend(list(map(str, Path(".").glob("ANTECHAMBER*"))))
                _ = list(map(lambda fname: fs.safe_rm(fname), files_to_remove))
                # gen frcmod
                # TODO(CJ): add some kind of check that this all actually runs correctly w/o errors
                run(
                    f"{Config.Amber.AmberHome}/bin/parmchk2 -i {prepin} -f prepi -o {frcmod}"
                )
                # record
                result.append((prepin, frcmod))
            return result

    def _combine_parm(
        self,
        lig_parms,
        prm_out_path="",
        o_dir="",
        ifsavepdb=0,
        ifsolve=1,
        box_type=None,
        box_size=3,  # Config.Amber.box_size,
        igb=None,
        if_prm_only=0,
    ):
        """
        combine different parmeter files and make finally inpcrd and prmtop
        -------
        structure: pdb
        ligands: prepi, frcmod
        metalcenters, artificial residues: TODO
        """
        if box_type == None:
            box_type = Config.Amber.box_type
        ll = list()
        leap_path = self.cache_path + "/leap.in"
        sol_path = self.path_name + "_ff.pdb"
        ll.extend(
            [
                "source leaprc.protein.ff14SB",
                "source leaprc.gaff",
                "source leaprc.water.tip3p",
            ]
        )
        # ligands
        for prepi, frcmod in lig_parms:
            ll.extend([f"loadAmberParams {frcmod}", f"loadAmberPrep {prepi}"])
        ll.append("a = loadpdb {self.path}")
        # igb Radii
        if igb != None:
            radii = radii_map[str(igb)]
            ll.append(f"set default PBRadii {radii}")
        of.write("center a" + line_feed)
        # solvation
        if ifsolve:
            of.write("addions a Na+ 0" + line_feed)
            of.write("addions a Cl- 0" + line_feed)
            if box_type == "box":
                of.write("solvatebox a TIP3PBOX " + box_size + line_feed)
            if box_type == "oct":
                of.write("solvateOct a TIP3PBOX " + box_size + line_feed)
            if box_type != "box" and box_type != "oct":
                raise Exception(
                    "PDB._combine_parm().box_type: Only support box and oct now!"
                )
        # save
        if prm_out_path == "":
            if o_dir == "":
                ll.append(
                    f"saveamberparm a {self.path_name}.prmtop {self.path_name}.inpcrd"
                )
                self.prmtop_path = self.path_name + ".prmtop"
                self.inpcrd_path = self.path_name + ".inpcrd"
            else:
                ll.append(
                    f"saveamberparm a {o_dir}/{self.name}.prmtop {o_dir}/{self.name}.inpcrd"
                )
                self.prmtop_path = o_dir + self.name + ".prmtop"
                self.inpcrd_path = o_dir + self.name + ".inpcrd"
        else:
            if o_dir == "":
                if if_prm_only:
                    mkdir("./tmp")
                    ll.append(f"saveamberparm a {prm_out_path} ./tmp/tmp.inpcrd")
                    self.prmtop_path = prm_out_path
                    self.inpcrd_path = None
                else:
                    ll.append(f"saveamberparm a {prm_out_path} {self.path_name}.inpcrd")
                    self.prmtop_path = prm_out_path
                    self.inpcrd_path = self.path_name + ".inpcrd"
            else:
                if if_prm_only:
                    mkdir("./tmp")
                    ll.append(f"saveamberparm a {prm_out_path} ./tmp/tmp.inpcrd")
                    self.prmtop_path = prm_out_path
                    self.inpcrd_path = None
                else:
                    of.write(
                        f"saveamberparm a {prm_out_path} {o_dir}/{self.name}.inpcrd"
                    )
                    self.prmtop_path = prm_out_path
                    self.inpcrd_path = o_dir + self.name + ".inpcrd"

        if ifsavepdb:
            of.write("savepdb a " + sol_path + line_feed)
        of.write("quit" + line_feed)
        # TODO(CJ): use EnvironmentManager() to run
        os.system("tleap -s -f " + leap_path + " > " + leap_path[:-2] + "out")

        return self.prmtop_path, self.inpcrd_path


#    def PDBMD(self, tag="", o_dir="", engine="Amber_GPU", equi_cpu=0):
#        """
#        Use self.prmtop_path and self.inpcrd_path to initilize a MD simulation.
#        The default MD configuration settings are assigned by class Config.Amber.
#        * User can also set MD configuration for the current object by assigning values in self.conf_xxxx.
#        * e.g.: self.conf_heat['nstlim'] = 50000
#        --------------
#        o_dir   : Write files in o_dir (current self.dir/MD by default).
#        tag     : tag the name of the MD folder
#        engine  : MD engine (cpu/gpu)
#        equi_cpu: if use cpu for equi step
#        Return the nc path of the prod step and store in self.nc
#        """
#        # make folder
#        if o_dir == "":
#            o_dir = self.dir + "/MD" + tag
#        mkdir(o_dir)
#
#        # express engine (pirority: AmberEXE_GPU/AmberEXE_CPU - AmberHome/bin/xxx)
#        PC_cmd, engine_path = Config.Amber.get_Amber_engine(engine=engine)
#        # express cpu engine if equi_cpu
#        if equi_cpu:
#            if Config.Amber.AmberEXE_CPU == None:
#                cpu_engine_path = Config.Amber.AmberHome + "/bin/sander.MPI"
#            else:
#                cpu_engine_path = Config.Amber.AmberEXE_CPU
#
#        # build input file (use self.MD_conf_xxxx)
#        min_path = self._build_MD_min(o_dir)
#        heat_path = self._build_MD_heat(o_dir)
#        equi_path = self._build_MD_equi(o_dir)
#        prod_path = self._build_MD_prod(o_dir)
#
#        # run sander
#        if Config.debug >= 1:
#            print(
#                "running: "
#                + PC_cmd
#                + " "
#                + engine_path
#                + " -O -i "
#                + min_path
#                + " -o "
#                + o_dir
#                + "/min.out -p "
#                + self.prmtop_path
#                + " -c "
#                + self.inpcrd_path
#                + " -r "
#                + o_dir
#                + "/min.rst -ref "
#                + self.inpcrd_path
#            )
#        os.system(
#            PC_cmd
#            + " "
#            + engine_path
#            + " -O -i "
#            + min_path
#            + " -o "
#            + o_dir
#            + "/min.out -p "
#            + self.prmtop_path
#            + " -c "
#            + self.inpcrd_path
#            + " -r "
#            + o_dir
#            + "/min.rst -ref "
#            + self.inpcrd_path
#        )
#        if Config.debug >= 1:
#            print(
#                "running: "
#                + PC_cmd
#                + " "
#                + engine_path
#                + " -O -i "
#                + heat_path
#                + " -o "
#                + o_dir
#                + "/heat.out -p "
#                + self.prmtop_path
#                + " -c "
#                + o_dir
#                + "/min.rst -ref "
#                + o_dir
#                + "/min.rst -r "
#                + o_dir
#                + "/heat.rst"
#            )
#        os.system(
#            PC_cmd
#            + " "
#            + engine_path
#            + " -O -i "
#            + heat_path
#            + " -o "
#            + o_dir
#            + "/heat.out -p "
#            + self.prmtop_path
#            + " -c "
#            + o_dir
#            + "/min.rst -ref "
#            + o_dir
#            + "/min.rst -r "
#            + o_dir
#            + "/heat.rst"
#        )
#
#        # gpu debug for equi
#        if equi_cpu:
#            # use Config.PC_cmd and cpu_engine_path
#            if Config.debug >= 1:
#                print(
#                    "running: "
#                    + Config.PC_cmd
#                    + " "
#                    + cpu_engine_path
#                    + " -O -i "
#                    + equi_path
#                    + " -o "
#                    + o_dir
#                    + "/equi.out -p "
#                    + self.prmtop_path
#                    + " -c "
#                    + o_dir
#                    + "/heat.rst -ref "
#                    + o_dir
#                    + "/heat.rst -r "
#                    + o_dir
#                    + "/equi.rst -x "
#                    + o_dir
#                    + "/equi.nc"
#                )
#            os.system(
#                Config.PC_cmd
#                + " "
#                + cpu_engine_path
#                + " -O -i "
#                + equi_path
#                + " -o "
#                + o_dir
#                + "/equi.out -p "
#                + self.prmtop_path
#                + " -c "
#                + o_dir
#                + "/heat.rst -ref "
#                + o_dir
#                + "/heat.rst -r "
#                + o_dir
#                + "/equi.rst -x "
#                + o_dir
#                + "/equi.nc"
#            )
#        else:
#            if Config.debug >= 1:
#                print(
#                    "running: "
#                    + PC_cmd
#                    + " "
#                    + engine_path
#                    + " -O -i "
#                    + equi_path
#                    + " -o "
#                    + o_dir
#                    + "/equi.out -p "
#                    + self.prmtop_path
#                    + " -c "
#                    + o_dir
#                    + "/heat.rst -ref "
#                    + o_dir
#                    + "/heat.rst -r "
#                    + o_dir
#                    + "/equi.rst -x "
#                    + o_dir
#                    + "/equi.nc"
#                )
#            os.system(
#                PC_cmd
#                + " "
#                + engine_path
#                + " -O -i "
#                + equi_path
#                + " -o "
#                + o_dir
#                + "/equi.out -p "
#                + self.prmtop_path
#                + " -c "
#                + o_dir
#                + "/heat.rst -ref "
#                + o_dir
#                + "/heat.rst -r "
#                + o_dir
#                + "/equi.rst -x "
#                + o_dir
#                + "/equi.nc"
#            )
#
#        if Config.debug >= 1:
#            print(
#                "running: "
#                + PC_cmd
#                + " "
#                + engine_path
#                + " -O -i "
#                + prod_path
#                + " -o "
#                + o_dir
#                + "/prod.out -p "
#                + self.prmtop_path
#                + " -c "
#                + o_dir
#                + "/equi.rst -ref "
#                + o_dir
#                + "/equi.rst -r "
#                + o_dir
#                + "/prod.rst -x "
#                + o_dir
#                + "/prod.nc"
#            )
#        os.system(
#            PC_cmd
#            + " "
#            + engine_path
#            + " -O -i "
#            + prod_path
#            + " -o "
#            + o_dir
#            + "/prod.out -p "
#            + self.prmtop_path
#            + " -c "
#            + o_dir
#            + "/equi.rst -ref "
#            + o_dir
#            + "/equi.rst -r "
#            + o_dir
#            + "/prod.rst -x "
#            + o_dir
#            + "/prod.nc"
#        )
#
#        self.nc = o_dir + "/prod.nc"
#        return o_dir + "/prod.nc"
#
#    def _build_MD_min(self, o_dir):
#        """
#        Build configuration file for a minimization job
#        See default value Config.Amber.conf_min
#        """
#        # path
#        o_path = o_dir + "/min.in"
#        # maxcyc related
#        maxcyc = self.conf_min["maxcyc"]
#        if self.conf_min["ncyc"] == "0.5maxcyc":
#            ncyc = str(int(0.5 * maxcyc))
#        if self.conf_min["ntpr"] == "0.01maxcyc":
#            ntpr = str(int(0.01 * maxcyc))
#        maxcyc = str(maxcyc)
#        # restrain related
#        if self.conf_min["ntr"] == "1":
#            ntr_line = (
#                "  ntr   = "
#                + self.conf_min["ntr"]
#                + ",	 restraint_wt = "
#                + self.conf_min["restraint_wt"]
#                + ", restraintmask = "
#                + self.conf_min["restraintmask"]
#                + ","
#                + line_feed
#            )
#        else:
#            ntr_line = ""
#        if self.conf_min["nmropt_rest"] == "1":
#            self.conf_min["nmropt"] = "1"
#            nmropt_line = "  nmropt= " + self.conf_min["nmropt"] + "," + line_feed
#            DISANG_tail = (
#                """ &wt
#    pe='END'
#
#    SANG= """
#                + self.conf_min["DISANG"]
#                + line_feed
#            )
#            self._build_MD_rs(step="min", o_path=self.conf_min["DISANG"])
#        else:
#            nmropt_line = ""
#            DISANG_tail = ""
#        # text
#        conf_str = (
#            """Minimize
#    trl
#    in  = 1,  ntx   = 1,  irest = 0,
#    c   = """
#            + self.conf_min["ntc"]
#            + """,    ntf = """
#            + self.conf_min["ntf"]
#            + """,
#    t   = """
#            + self.conf_min["cut"]
#            + """,
#    xcyc= """
#            + maxcyc
#            + """, ncyc  = """
#            + ncyc
#            + """,
#    pr  = """
#            + ntpr
#            + """, ntwx  = 0,
#
#            + ntr_line
#            + nmropt_line
#            + """ /
#
#            + DISANG_tail
#        )
#        # write
#        with open(o_path, "w") as of:
#            of.write(conf_str)
#        return o_path
#
#    def _build_MD_heat(self, o_dir):
#        """
#        Build configuration file for a heat job
#        See default value Config.Amber.conf_heat
#        """
#        # path
#        o_path = o_dir + "/heat.in"
#
#        # nstlim related
#        nstlim = self.conf_heat["nstlim"]
#        if self.conf_heat["A_istep2"] == "0.9nstlim":
#            A_istep2 = str(int(nstlim * 0.9))
#        if self.conf_heat["B_istep1"] == "A_istep2+1":
#            B_istep1 = str(int(A_istep2) + 1)
#        if self.conf_heat["ntpr"] == "0.01nstlim":
#            ntpr = str(int(nstlim * 0.01))
#        if self.conf_heat["ntwx"] == "nstlim":
#            ntwx = str(nstlim)
#        nstlim = str(nstlim)
#        # restrain related
#        if self.conf_heat["ntr"] == "1":
#            ntr_line = (
#                "  ntr   = "
#                + self.conf_heat["ntr"]
#                + ", restraint_wt = "
#                + self.conf_heat["restraint_wt"]
#                + ", restraintmask = "
#                + self.conf_heat["restraintmask"]
#                + ","
#                + line_feed
#            )
#        else:
#            ntr_line = ""
#        if self.conf_heat["nmropt_rest"] == "1":
#            DISANG_tail = """  DISANG=""" + self.conf_heat["DISANG"] + line_feed
#            self._build_MD_rs(step="heat", o_path=self.conf_heat["DISANG"])
#        else:
#            DISANG_tail = ""
#        conf_str = (
#            '''Heat
#    trl
#    in  = 0,  ntx = 1, irest = 0,
#    c   = """
#            + self.conf_heat["ntc"]
#            + """, ntf = """
#            + self.conf_heat["ntf"]
#            + """,
#    t   = """
#            + self.conf_heat["cut"]
#            + """,
#    tlim= """
#            + nstlim
#            + """, dt= """
#            + self.conf_heat["dt"]
#            + """,
#    mpi = """
#            + self.conf_heat["tempi"]
#            + """,  temp0="""
#            + self.conf_heat["temp0"]
#            + """,
#    pr  = """
#            + ntpr
#            + """,  ntwx="""
#            + ntwx
#            + """,
#    t   = """
#            + self.conf_heat["ntt"]
#            + """, gamma_ln = """
#            + self.conf_heat["gamma_ln"]
#            + """,
#    b   = 1,  ntp = 0,
#    rap = """
#            + self.conf_heat["iwarp"]
#            + """,
#    ropt= 1,
#        = -1,
#
#            + ntr_line
#            + """ /
#
#    pe  = 'TEMP0',
#    tep1= 0, istep2="""
#            + A_istep2
#            + """,
#    lue1= """
#            + self.conf_heat["tempi"]
#            + """, value2="""
#            + self.conf_heat["temp0"]
#            + """,
#
#
#    pe  = 'TEMP0',
#    tep1= """
#            + B_istep1
#            + """, istep2="""
#            + nstlim
#            + """,
#    lue1= """
#            + self.conf_heat["temp0"]
#            + """, value2="""
#            + self.conf_heat["temp0"]
#            + """,
#
#
#    pe  = 'END',
#
#
#            + DISANG_tail
#        )
#        # write
#        with open(o_path, "w") as of:
#            of.write(conf_str)
#        return o_path
#
#    def _build_MD_equi(self, o_dir):
#        """
#        Build configuration file for the equilibration step
#        See default value Config.Amber.conf_equi
#        default ntwx -> 10ps
#        """
#        # path
#        o_path = o_dir + "/equi.in"
#        #
#        #        # nstlim related
#        #        nstlim = self.conf_equi["nstlim"]
#        #        if self.conf_equi["ntpr"] == "0.002nstlim":
#        #            ntpr = str(int(nstlim * 0.002))
#        #        nstlim = str(nstlim)
#        #        # restrain related
#        #        if self.conf_equi["ntr"] == "1":
#        #            ntr_line = (
#        #                "  ntr   = "
#        #                + self.conf_equi["ntr"]
#        #                + ", restraint_wt = "
#        #                + self.conf_equi["restraint_wt"]
#        #                + ", restraintmask = "
#        #                + self.conf_equi["restraintmask"]
#        #                + ","
#        #                + line_feed
#        #            )
#        #        else:
#        #            ntr_line = ''
#        #        if self.conf_equi['nmropt_rest'] == '1':
#        #            self.conf_equi['nmropt'] = '1'
#        #            nmropt_line = '  nmropt= '+self.conf_equi['nmropt']+','+line_feed
#        #            DISANG_tail = ''' &wt
#        #  type='END'
#        # /
#        #  DISANG= '''+self.conf_equi['DISANG']+line_feed
#        #            self._build_MD_rs(step='equi',o_path=self.conf_equi['DISANG'])
#        #        else:
#        #            nmropt_line = ''
#        #            DISANG_tail = ''
#        #
#        #        conf_str = (
#        #            """Equilibration:constant pressure
#        # &cntrl
#        #  imin  = 0,  ntx = """
#        #            + self.conf_equi["ntx"]
#        #            + """,  irest = """
#        #            + self.conf_equi["irest"]
#        #            + """,
#        #  ntf   = """
#        #            + self.conf_equi["ntf"]
#        #            + """,  ntc = """
#        #            + self.conf_equi["ntc"]
#        #            + """,
#        #  nstlim= """
#        #            + nstlim
#        #            + """, dt= """
#        #            + self.conf_equi["dt"]
#        #            + """,
#        #  cut   = """
#        #            + self.conf_equi["cut"]
#        #            + """,
#        #  temp0 = """
#        #            + self.conf_equi["temp0"]
#        #            + """,
#        #  ntpr  = """
#        #            + ntpr
#        #            + """, ntwx = """
#        #            + self.conf_equi["ntwx"]
#        #            + """,
#        #  ntt   = """
#        #            + self.conf_equi["ntt"]
#        #            + """, gamma_ln = """
#        #            + self.conf_equi["gamma_ln"]
#        #            + """,
#        #  ntb   = 2,  ntp = 1,
#        #  iwrap = """
#        #            + self.conf_equi["iwarp"]
#        #            + """,
#        #  ig    = -1,
#        # """+ntr_line+nmropt_line+''' /
#        #'''+DISANG_tail
#        #        #write
#        with open(o_path, "w") as of:
#            of.write(conf_str)
#        return o_path
#
#    def _build_MD_prod(self, o_dir):
#        """
#        Build configuration file for the production step
#        See default value Config.Amber.conf_prod
#        default ntwx -> 10ps
#        """
#        # path
#        o_path = o_dir + "/prod.in"
#
#        #        # nstlim related
#        #        nstlim = self.conf_prod["nstlim"]
#        #        if self.conf_prod["ntpr"] == "0.001nstlim":
#        #            ntpr = str(int(nstlim * 0.001))
#        #        nstlim = str(nstlim)
#        #        # restrain related
#        #        if self.conf_prod["ntr"] == "1":
#        #            ntr_line = (
#        #                "  ntr   = "
#        #                + self.conf_prod["ntr"]
#        #                + ", restraint_wt = "
#        #                + self.conf_prod["restraint_wt"]
#        #                + ", restraintmask = "
#        #                + self.conf_prod["restraintmask"]
#        #                + ","
#        #                + line_feed
#        #            )
#        #        else:
#        #            ntr_line = ''
#        #        if self.conf_prod['nmropt_rest'] == '1':
#        #            self.conf_prod['nmropt'] = '1'
#        #            nmropt_line = '  nmropt= '+self.conf_prod['nmropt']+','+line_feed
#        #            DISANG_tail = ''' &wt
#        #  type='END'
#        # /
#        #  DISANG= '''+self.conf_prod['DISANG']+line_feed
#        #            self._build_MD_rs(step='prod',o_path=self.conf_prod['DISANG'])
#        #        else:
#        #            nmropt_line = ''
#        #            DISANG_tail = ''
#        #
#        #        conf_str = (
#        #            """Production: constant pressure
#        # &cntrl
#        #  imin  = 0, ntx = """
#        #            + self.conf_prod["ntx"]
#        #            + """, irest = """
#        #            + self.conf_prod["irest"]
#        #            + """,
#        #  ntf   = """
#        #            + self.conf_prod["ntf"]
#        #            + """,  ntc = """
#        #            + self.conf_prod["ntc"]
#        #            + """,
#        #  nstlim= """
#        #            + nstlim
#        #            + """, dt= """
#        #            + self.conf_prod["dt"]
#        #            + """,
#        #  cut   = """
#        #            + self.conf_prod["cut"]
#        #            + """,
#        #  temp0 = """
#        #            + self.conf_prod["temp0"]
#        #            + """,
#        #  ntpr  = """
#        #            + ntpr
#        #            + """, ntwx = """
#        #            + self.conf_prod["ntwx"]
#        #            + """,
#        #  ntt   = """
#        #            + self.conf_prod["ntt"]
#        #            + """, gamma_ln = """
#        #            + self.conf_prod["gamma_ln"]
#        #            + """,
#        #  ntb   = 2,  ntp = 1,
#        #  iwrap = """
#        #            + self.conf_prod["iwarp"]
#        #            + """,
#        #  ig    = -1,
#        # """+ntr_line+nmropt_line+''' /
#        #'''+DISANG_tail
#        # write
#        with open(o_path, "w") as of:
#            of.write(conf_str)
#        return o_path
#
#    def _build_MD_rs(self, step, o_path):
#        """
#        Generate a file for DISANG restraint. Get parameters from self.conf_step.
#        """
#        rs_str = ""
#        rs_data_step = self.__dict__["conf_" + step]["rs_constraints"]
#        for rest_data in rs_data_step:
#            rs_str = (
#                rs_str
#                + """  &rst
#    at=  """
#                + ",".join(rest_data["iat"])
#                + """, r1= """
#                + rest_data["r1"]
#                + """, r2= """
#                + rest_data["r2"]
#                + """, r3= """
#                + rest_data["r3"]
#                + """, r4= """
#                + rest_data["r4"]
#                + """,
#    k2="""
#                + rest_data["rk2"]
#                + """, rk3="""
#                + rest_data["rk3"]
#                + """, ir6="""
#                + rest_data["ir6"]
#                + """, ialtd="""
#                + rest_data["ialtd"]
#                + """,
#    nd
#
#            )
#        # write
#        with open(o_path, "w") as of:
#            of.write(rs_str)
#        return o_path
#
#    def reset_MD_conf(self):
#        """
#        reset MD configuration of current object to default.
#        """
#        self.conf_min = Config.Amber.conf_min
#        self.conf_heat = Config.Amber.conf_heat
#        self.conf_equi = Config.Amber.conf_equi
#        self.conf_prod = Config.Amber.conf_prod
#
#    def show_MD_conf(self):
#        """
#        Show MD configuration of current object.
#        """
#        print("Min :     " + repr(self.conf_min))
#        print("Heat:     " + repr(self.conf_heat))
#        print("Equi:     " + repr(self.conf_equi))
#        print("Prod:     " + repr(self.conf_prod))
#
#    def nc2mdcrd(
#        self, o_path="", point=None, start=1, end=-1, step=1, engine="cpptraj"
#    ):
#        """
#        convert self.nc to a mdcrd file to read and operate.(self.nc[:-2]+'.mdcrd' by default)
#        a easier way is to use pytraj directly.
#        ---------------
#        o_path: user assigned out path (self.nc[:-2]+'mdcrd' by default)
#        point:  sample point. use value from self.conf_prod['nstlim'] and self.conf_prod['ntwx'] to determine step size.
#        start:  start point
#        end:    end point
#        step:   step size
#        engine: pytraj or cpptraj (some package conflict may cause pytraj not available)
#        """
#        if self.nc == None:
#            raise Exception(
#                "No nc file found. Please assign self.nc or run PDBMD first"
#            )
#        else:
#            if o_path == "":
#                o_path = self.nc[:-2] + "mdcrd"
#            if end == -1:
#                end = "last"
#            if point != None:
#                all_p = int(self.conf_prod["nstlim"]) / int(self.conf_prod["ntwx"])
#                step = int(all_p / point)
#
#            if engine not in ["pytraj", "cpptraj"]:
#                raise Exception("engine: pytraj or cpptraj")
#
#            if engine == "pytraj":
#                pass
#
#            if engine == "cpptraj":
#                cpp_in_path = self.cache_path + "/cpptraj_nc2mdcrd.in"
#                cpp_out_path = self.cache_path + "/cpptraj_nc2mdcrd.out"
#                with open(cpp_in_path, "w") as of:
#                    of.write("parm " + self.prmtop_path + line_feed)
#                    of.write(
#                        "trajin "
#                        + self.nc
#                        + " "
#                        + str(start)
#                        + " "
#                        + end
#                        + " "
#                        + str(step)
#                        + line_feed
#                    )
#                    of.write("trajout " + o_path + line_feed)
#                    of.write("run" + line_feed)
#                    of.write("quit" + line_feed)
#                os.system("cpptraj -i " + cpp_in_path + " > " + cpp_out_path)
#
#        self.mdcrd = o_path
#        return o_path

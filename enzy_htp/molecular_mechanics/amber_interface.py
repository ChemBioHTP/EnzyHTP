# TODO documentation
import shutil
from ..core.logger import _LOGGER

class AmberInterface:
    def __init__(self, config=None):
        self._confg = config

    def minimize_structure(self, pdb : str, mode: str = 'CPU' ):

    #def PDBMin(self, cycle=2000, engine="Amber_GPU"):
        """
        Run a minization use self.prmtop and self.inpcrd and setting form class Config.
        --------------------------------------------------
        Save changed PDB to self.path (containing water and ions)
        Mainly used for remove bad contact from PDB2PDBwLeap.
        """

        out4_PDB_path = self.path_name + "_min.pdb"

        # make sander input
        min_dir = self.cache_path + "/PDBMin"
        minin_path = min_dir + "/min.in"
        minout_path = min_dir + "/min.out"
        minrst_path = (min_dir + "/min.ncrst"
                      )  # change to ncrst seeking for solution of rst error
        mkdir(min_dir)
        min_input = open(minin_path, "w")
        min_input.write("Minimize" + line_feed)
        min_input.write(" &cntrl" + line_feed)
        min_input.write("  imin=1," + line_feed)
        min_input.write("  ntx=1," + line_feed)
        min_input.write("  irest=0," + line_feed)
        min_input.write("  maxcyc=" + str(cycle) + "," + line_feed)
        min_input.write("  ncyc=" + str(int(0.5 * cycle)) + "," + line_feed)
        min_input.write("  ntpr=" + str(int(0.2 * cycle)) + "," + line_feed)
        min_input.write("  ntwx=0," + line_feed)
        min_input.write("  cut=8.0," + line_feed)
        min_input.write(" /" + line_feed)
        min_input.close()

        # express engine
        PC_cmd, engine_path = Config.Amber.get_Amber_engine(engine=engine)
        engine = self._config.amber_confg()[ 'CPU_ENGINE' if mode == 'CPU' else 'GPU_ENGINE' ]
        self._config.env_manager().run_command( engine, f"-O -i {minin_path} -o {minout_path} -p {self.prmtop_path} -c {self.inpcrd_path} -r {minrst_path}")
        self._config.env_manager().run_command('ambpdp', f"-p {self.prmtop_path} -c {minrst_path} > {out4_PDB_path}")

        os.system("mv " + self.prmtop_path + " " + self.inpcrd_path + " " +
                  min_dir) # TODO CJ: not sure if this actually works...

        self.path = out4_PDB_path
        self._update_name()
        self.prmtop_path = min_dir + "/" + self.prmtop_path
        self.inpcrd_path = min_dir + "/" + self.inpcrd_path


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

        ligands_pathNchrg = self.stru.build_ligands(lig_dir,
                                                    ifcharge=1,
                                                    ifunique=1)
        # metalcenters_path = self.stru.build_metalcenters(met_dir)
        # parm
        ligand_parm_paths = self._ligand_parm(ligands_pathNchrg,
                                              method=lig_method,
                                              renew=renew_lig)
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

    def _ligand_parm(self, paths, method="AM1BCC", renew=0):
        """
        Turn ligands to prepi (w/net charge), parameterize with parmchk
        return [(perpi_1, frcmod_1), ...]
        -----------
        method  : method use for ligand charge. Only support AM1BCC now.
        renew   : 0:(default) use old parm files if exist. 1: renew parm files everytime
        TODO check if the ligand is having correct name. (isolate the renaming function and also use in the class structure)
        * WARN: The parm file for ligand will always be like xxx/ligand_1.frcmod. Remember to enable renew when different object is sharing a same path.
        * BUG: Antechamber has a bug that if current dir has temp files from previous antechamber run (ANTECHAMBER_AC.AC, etc.) sqm will fail. Now remove them everytime.
        """
        parm_paths = []
        self.prepi_path = {}

        for lig_pdb, net_charge in paths:
            if method == "AM1BCC":
                out_prepi = lig_pdb[:-3] + "prepin"
                out_frcmod = lig_pdb[:-3] + "frcmod"
                with open(lig_pdb) as f:
                    for line in f:
                        pdbl = PDB_line(line)
                        if pdbl.line_type == "ATOM" or pdbl.line_type == "HETATM":
                            lig_name = pdbl.resi_name
                # if renew
                if (os.path.isfile(out_prepi) and os.path.isfile(out_frcmod) and
                        not renew):
                    if Config.debug >= 1:
                        print("Parm files exist: " + out_prepi + " " +
                              out_frcmod)
                        print("Using old parm files.")
                else:
                    # gen prepi (net charge and correct protonation state is important)
                    if Config.debug >= 1:
                        print("running: " + Config.Amber.AmberHome +
                              "/bin/antechamber -i " + lig_pdb +
                              " -fi pdb -o " + out_prepi +
                              " -fo prepi -c bcc -s 0 -nc " + str(net_charge))
                    run(
                        Config.Amber.AmberHome + "/bin/antechamber -i " +
                        lig_pdb + " -fi pdb -o " + out_prepi +
                        " -fo prepi -c bcc -s 0 -nc " + str(net_charge),
                        check=True,
                        text=True,
                        shell=True,
                        capture_output=True,
                    )
                    if Config.debug <= 1:
                        os.system(
                            "rm ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.pdb sqm.in sqm.out"
                        )
                    # gen frcmod
                    if Config.debug >= 1:
                        print("running: " + Config.Amber.AmberHome +
                              "/bin/parmchk2 -i " + out_prepi +
                              " -f prepi -o " + out_frcmod)
                    run(
                        Config.Amber.AmberHome + "/bin/parmchk2 -i " +
                        out_prepi + " -f prepi -o " + out_frcmod,
                        check=True,
                        text=True,
                        shell=True,
                        capture_output=True,
                    )
                # record
                parm_paths.append((out_prepi, out_frcmod))
                self.prepi_path[lig_name] = out_prepi

        return parm_paths

        return self.path

        pass


    def _combine_parm(
        self,
        lig_parms,
        prm_out_path="",
        o_dir="",
        ifsavepdb=0,
        ifsolve=1,
        box_type=None,
        box_size=3,  #Config.Amber.box_size,
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

        leap_path = self.cache_path + "/leap.in"
        sol_path = self.path_name + "_ff.pdb"
        with open(leap_path, "w") as of:
            of.write("source leaprc.protein.ff14SB" + line_feed)
            of.write("source leaprc.gaff" + line_feed)
            of.write("source leaprc.water.tip3p" + line_feed)
            # ligands
            for prepi, frcmod in lig_parms:
                of.write("loadAmberParams " + frcmod + line_feed)
                of.write("loadAmberPrep " + prepi + line_feed)
            of.write("a = loadpdb " + self.path + line_feed)
            # igb Radii
            if igb != None:
                radii = radii_map[str(igb)]
                of.write("set default PBRadii " + radii + line_feed)
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
                    of.write("saveamberparm a " + self.path_name + ".prmtop " +
                             self.path_name + ".inpcrd" + line_feed)
                    self.prmtop_path = self.path_name + ".prmtop"
                    self.inpcrd_path = self.path_name + ".inpcrd"
                else:
                    of.write("saveamberparm a " + o_dir + self.name +
                             ".prmtop " + o_dir + self.name + ".inpcrd" +
                             line_feed)
                    self.prmtop_path = o_dir + self.name + ".prmtop"
                    self.inpcrd_path = o_dir + self.name + ".inpcrd"
            else:
                if o_dir == "":
                    if if_prm_only:
                        mkdir("./tmp")
                        of.write("saveamberparm a " + prm_out_path +
                                 " ./tmp/tmp.inpcrd" + line_feed)
                        self.prmtop_path = prm_out_path
                        self.inpcrd_path = None
                    else:
                        of.write("saveamberparm a " + prm_out_path + " " +
                                 self.path_name + ".inpcrd" + line_feed)
                        self.prmtop_path = prm_out_path
                        self.inpcrd_path = self.path_name + ".inpcrd"
                else:
                    if if_prm_only:
                        mkdir("./tmp")
                        of.write("saveamberparm a " + prm_out_path +
                                 " ./tmp/tmp.inpcrd" + line_feed)
                        self.prmtop_path = prm_out_path
                        self.inpcrd_path = None
                    else:
                        of.write("saveamberparm a " + prm_out_path + " " +
                                 o_dir + self.name + ".inpcrd" + line_feed)
                        self.prmtop_path = prm_out_path
                        self.inpcrd_path = o_dir + self.name + ".inpcrd"

            if ifsavepdb:
                of.write("savepdb a " + sol_path + line_feed)
            of.write("quit" + line_feed)

        os.system("tleap -s -f " + leap_path + " > " + leap_path[:-2] + "out")

        return self.prmtop_path, self.inpcrd_path


    def PDBMD(self, tag="", o_dir="", engine="Amber_GPU", equi_cpu=0):
        """
        Use self.prmtop_path and self.inpcrd_path to initilize a MD simulation.
        The default MD configuration settings are assigned by class Config.Amber.
        * User can also set MD configuration for the current object by assigning values in self.conf_xxxx.
        * e.g.: self.conf_heat['nstlim'] = 50000
        --------------
        o_dir   : Write files in o_dir (current self.dir/MD by default).
        tag     : tag the name of the MD folder
        engine  : MD engine (cpu/gpu)
        equi_cpu: if use cpu for equi step
        Return the nc path of the prod step and store in self.nc
        """
        # make folder
        if o_dir == "":
            o_dir = self.dir + "/MD" + tag
        mkdir(o_dir)

        # express engine (pirority: AmberEXE_GPU/AmberEXE_CPU - AmberHome/bin/xxx)
        PC_cmd, engine_path = Config.Amber.get_Amber_engine(engine=engine)
        # express cpu engine if equi_cpu
        if equi_cpu:
            if Config.Amber.AmberEXE_CPU == None:
                cpu_engine_path = Config.Amber.AmberHome + "/bin/sander.MPI"
            else:
                cpu_engine_path = Config.Amber.AmberEXE_CPU

        # build input file (use self.MD_conf_xxxx)
        min_path = self._build_MD_min(o_dir)
        heat_path = self._build_MD_heat(o_dir)
        equi_path = self._build_MD_equi(o_dir)
        prod_path = self._build_MD_prod(o_dir)

        # run sander
        if Config.debug >= 1:
            print("running: " + PC_cmd + " " + engine_path + " -O -i " +
                  min_path + " -o " + o_dir + "/min.out -p " +
                  self.prmtop_path + " -c " + self.inpcrd_path + " -r " +
                  o_dir + "/min.rst -ref " + self.inpcrd_path)
        os.system(PC_cmd + " " + engine_path + " -O -i " + min_path + " -o " +
                  o_dir + "/min.out -p " + self.prmtop_path + " -c " +
                  self.inpcrd_path + " -r " + o_dir + "/min.rst -ref " +
                  self.inpcrd_path)
        if Config.debug >= 1:
            print("running: " + PC_cmd + " " + engine_path + " -O -i " +
                  heat_path + " -o " + o_dir + "/heat.out -p " +
                  self.prmtop_path + " -c " + o_dir + "/min.rst -ref " + o_dir +
                  "/min.rst -r " + o_dir + "/heat.rst")
        os.system(PC_cmd + " " + engine_path + " -O -i " + heat_path + " -o " +
                  o_dir + "/heat.out -p " + self.prmtop_path + " -c " + o_dir +
                  "/min.rst -ref " + o_dir + "/min.rst -r " + o_dir +
                  "/heat.rst")

        # gpu debug for equi
        if equi_cpu:
            # use Config.PC_cmd and cpu_engine_path
            if Config.debug >= 1:
                print("running: " + Config.PC_cmd + " " + cpu_engine_path +
                      " -O -i " + equi_path + " -o " + o_dir + "/equi.out -p " +
                      self.prmtop_path + " -c " + o_dir + "/heat.rst -ref " +
                      o_dir + "/heat.rst -r " + o_dir + "/equi.rst -x " +
                      o_dir + "/equi.nc")
            os.system(Config.PC_cmd + " " + cpu_engine_path + " -O -i " +
                      equi_path + " -o " + o_dir + "/equi.out -p " +
                      self.prmtop_path + " -c " + o_dir + "/heat.rst -ref " +
                      o_dir + "/heat.rst -r " + o_dir + "/equi.rst -x " +
                      o_dir + "/equi.nc")
        else:
            if Config.debug >= 1:
                print("running: " + PC_cmd + " " + engine_path + " -O -i " +
                      equi_path + " -o " + o_dir + "/equi.out -p " +
                      self.prmtop_path + " -c " + o_dir + "/heat.rst -ref " +
                      o_dir + "/heat.rst -r " + o_dir + "/equi.rst -x " +
                      o_dir + "/equi.nc")
            os.system(PC_cmd + " " + engine_path + " -O -i " + equi_path +
                      " -o " + o_dir + "/equi.out -p " + self.prmtop_path +
                      " -c " + o_dir + "/heat.rst -ref " + o_dir +
                      "/heat.rst -r " + o_dir + "/equi.rst -x " + o_dir +
                      "/equi.nc")

        if Config.debug >= 1:
            print("running: " + PC_cmd + " " + engine_path + " -O -i " +
                  prod_path + " -o " + o_dir + "/prod.out -p " +
                  self.prmtop_path + " -c " + o_dir + "/equi.rst -ref " +
                  o_dir + "/equi.rst -r " + o_dir + "/prod.rst -x " + o_dir +
                  "/prod.nc")
        os.system(PC_cmd + " " + engine_path + " -O -i " + prod_path + " -o " +
                  o_dir + "/prod.out -p " + self.prmtop_path + " -c " + o_dir +
                  "/equi.rst -ref " + o_dir + "/equi.rst -r " + o_dir +
                  "/prod.rst -x " + o_dir + "/prod.nc")

        self.nc = o_dir + "/prod.nc"
        return o_dir + "/prod.nc"


    def _build_MD_min(self, o_dir):
        """
        Build configuration file for a minimization job
        See default value Config.Amber.conf_min
        """
        # path
        o_path = o_dir + "/min.in"
        # maxcyc related
        maxcyc = self.conf_min["maxcyc"]
        if self.conf_min["ncyc"] == "0.5maxcyc":
            ncyc = str(int(0.5 * maxcyc))
        if self.conf_min["ntpr"] == "0.01maxcyc":
            ntpr = str(int(0.01 * maxcyc))
        maxcyc = str(maxcyc)
        # restrain related
        if self.conf_min["ntr"] == "1":
            ntr_line = ("  ntr   = " + self.conf_min["ntr"] +
                        ",	 restraint_wt = " + self.conf_min["restraint_wt"] +
                        ", restraintmask = " + self.conf_min["restraintmask"] +
                        "," + line_feed)
        else:
            ntr_line = ''
        if self.conf_min['nmropt_rest'] == '1':
            self.conf_min['nmropt'] = '1'
            nmropt_line = '  nmropt= ' + self.conf_min[
                'nmropt'] + ',' + line_feed
            DISANG_tail = ''' &wt
  type='END'
 /
  DISANG= ''' + self.conf_min['DISANG'] + line_feed
            self._build_MD_rs(step='min', o_path=self.conf_min['DISANG'])
        else:
            nmropt_line = ''
            DISANG_tail = ''
        #text
        conf_str = '''Minimize
 &cntrl
  imin  = 1,  ntx   = 1,  irest = 0,
  ntc   = ''' + self.conf_min['ntc'] + ''',    ntf = ''' + self.conf_min[
            'ntf'] + ''',
  cut   = ''' + self.conf_min['cut'] + ''',
  maxcyc= ''' + maxcyc + ''', ncyc  = ''' + ncyc + ''',
  ntpr  = ''' + ntpr + ''', ntwx  = 0,
''' + ntr_line + nmropt_line + ''' /
''' + DISANG_tail
        #write
        with open(o_path, 'w') as of:
            of.write(conf_str)
        return o_path

    def _build_MD_heat(self, o_dir):
        """
        Build configuration file for a heat job
        See default value Config.Amber.conf_heat
        """
        # path
        o_path = o_dir + "/heat.in"

        # nstlim related
        nstlim = self.conf_heat["nstlim"]
        if self.conf_heat["A_istep2"] == "0.9nstlim":
            A_istep2 = str(int(nstlim * 0.9))
        if self.conf_heat["B_istep1"] == "A_istep2+1":
            B_istep1 = str(int(A_istep2) + 1)
        if self.conf_heat["ntpr"] == "0.01nstlim":
            ntpr = str(int(nstlim * 0.01))
        if self.conf_heat["ntwx"] == "nstlim":
            ntwx = str(nstlim)
        nstlim = str(nstlim)
        # restrain related
        if self.conf_heat["ntr"] == "1":
            ntr_line = ("  ntr   = " + self.conf_heat["ntr"] +
                        ", restraint_wt = " + self.conf_heat["restraint_wt"] +
                        ", restraintmask = " + self.conf_heat["restraintmask"] +
                        "," + line_feed)
        else:
            ntr_line = ''
        if self.conf_heat['nmropt_rest'] == '1':
            DISANG_tail = '''  DISANG=''' + self.conf_heat['DISANG'] + line_feed
            self._build_MD_rs(step='heat', o_path=self.conf_heat['DISANG'])
        else:
            DISANG_tail = ''
        conf_str = '''Heat
 &cntrl
  imin  = 0,  ntx = 1, irest = 0,
  ntc   = """
            + self.conf_heat["ntc"]
            + """, ntf = """
            + self.conf_heat["ntf"]
            + """,
  cut   = """
            + self.conf_heat["cut"]
            + """,
  nstlim= """
            + nstlim
            + """, dt= """
            + self.conf_heat["dt"]
            + """,
  tempi = """
            + self.conf_heat["tempi"]
            + """,  temp0="""
            + self.conf_heat["temp0"]
            + """,  
  ntpr  = """
            + ntpr
            + """,  ntwx="""
            + ntwx
            + """,
  ntt   = """
            + self.conf_heat["ntt"]
            + """, gamma_ln = """
            + self.conf_heat["gamma_ln"]
            + """,
  ntb   = 1,  ntp = 0,
  iwrap = """
            + self.conf_heat["iwarp"]
            + """,
  nmropt= 1,
  ig    = -1,
"""
            + ntr_line
            + """ /
 &wt
  type  = 'TEMP0',
  istep1= 0, istep2="""
            + A_istep2
            + """,
  value1= """
            + self.conf_heat["tempi"]
            + """, value2="""
            + self.conf_heat["temp0"]
            + """,
 /
 &wt
  type  = 'TEMP0',
  istep1= """
            + B_istep1
            + """, istep2="""
            + nstlim
            + """,
  value1= """
            + self.conf_heat["temp0"]
            + """, value2="""
            + self.conf_heat["temp0"]
            + """,
 /
 &wt
  type  = 'END',
 /
''' + DISANG_tail
        #write
        with open(o_path, 'w') as of:
            of.write(conf_str)
        return o_path

    def _build_MD_equi(self, o_dir):
        """
        Build configuration file for the equilibration step
        See default value Config.Amber.conf_equi
        default ntwx -> 10ps
        """
        # path
        o_path = o_dir + "/equi.in"
        #
        #        # nstlim related
        #        nstlim = self.conf_equi["nstlim"]
        #        if self.conf_equi["ntpr"] == "0.002nstlim":
        #            ntpr = str(int(nstlim * 0.002))
        #        nstlim = str(nstlim)
        #        # restrain related
        #        if self.conf_equi["ntr"] == "1":
        #            ntr_line = (
        #                "  ntr   = "
        #                + self.conf_equi["ntr"]
        #                + ", restraint_wt = "
        #                + self.conf_equi["restraint_wt"]
        #                + ", restraintmask = "
        #                + self.conf_equi["restraintmask"]
        #                + ","
        #                + line_feed
        #            )
        #        else:
        #            ntr_line = ''
        #        if self.conf_equi['nmropt_rest'] == '1':
        #            self.conf_equi['nmropt'] = '1'
        #            nmropt_line = '  nmropt= '+self.conf_equi['nmropt']+','+line_feed
        #            DISANG_tail = ''' &wt
        #  type='END'
        # /
        #  DISANG= '''+self.conf_equi['DISANG']+line_feed
        #            self._build_MD_rs(step='equi',o_path=self.conf_equi['DISANG'])
        #        else:
        #            nmropt_line = ''
        #            DISANG_tail = ''
        #
        #        conf_str = (
        #            """Equilibration:constant pressure
        # &cntrl
        #  imin  = 0,  ntx = """
        #            + self.conf_equi["ntx"]
        #            + """,  irest = """
        #            + self.conf_equi["irest"]
        #            + """,
        #  ntf   = """
        #            + self.conf_equi["ntf"]
        #            + """,  ntc = """
        #            + self.conf_equi["ntc"]
        #            + """,
        #  nstlim= """
        #            + nstlim
        #            + """, dt= """
        #            + self.conf_equi["dt"]
        #            + """,
        #  cut   = """
        #            + self.conf_equi["cut"]
        #            + """,
        #  temp0 = """
        #            + self.conf_equi["temp0"]
        #            + """,
        #  ntpr  = """
        #            + ntpr
        #            + """, ntwx = """
        #            + self.conf_equi["ntwx"]
        #            + """,
        #  ntt   = """
        #            + self.conf_equi["ntt"]
        #            + """, gamma_ln = """
        #            + self.conf_equi["gamma_ln"]
        #            + """,
        #  ntb   = 2,  ntp = 1,
        #  iwrap = """
        #            + self.conf_equi["iwarp"]
        #            + """,
        #  ig    = -1,
        #"""+ntr_line+nmropt_line+''' /
        #'''+DISANG_tail
        #        #write
        with open(o_path, 'w') as of:
            of.write(conf_str)
        return o_path

    def _build_MD_prod(self, o_dir):
        """
        Build configuration file for the production step
        See default value Config.Amber.conf_prod
        default ntwx -> 10ps
        """
        # path
        o_path = o_dir + "/prod.in"

        #        # nstlim related
        #        nstlim = self.conf_prod["nstlim"]
        #        if self.conf_prod["ntpr"] == "0.001nstlim":
        #            ntpr = str(int(nstlim * 0.001))
        #        nstlim = str(nstlim)
        #        # restrain related
        #        if self.conf_prod["ntr"] == "1":
        #            ntr_line = (
        #                "  ntr   = "
        #                + self.conf_prod["ntr"]
        #                + ", restraint_wt = "
        #                + self.conf_prod["restraint_wt"]
        #                + ", restraintmask = "
        #                + self.conf_prod["restraintmask"]
        #                + ","
        #                + line_feed
        #            )
        #        else:
        #            ntr_line = ''
        #        if self.conf_prod['nmropt_rest'] == '1':
        #            self.conf_prod['nmropt'] = '1'
        #            nmropt_line = '  nmropt= '+self.conf_prod['nmropt']+','+line_feed
        #            DISANG_tail = ''' &wt
        #  type='END'
        # /
        #  DISANG= '''+self.conf_prod['DISANG']+line_feed
        #            self._build_MD_rs(step='prod',o_path=self.conf_prod['DISANG'])
        #        else:
        #            nmropt_line = ''
        #            DISANG_tail = ''
        #
        #        conf_str = (
        #            """Production: constant pressure
        # &cntrl
        #  imin  = 0, ntx = """
        #            + self.conf_prod["ntx"]
        #            + """, irest = """
        #            + self.conf_prod["irest"]
        #            + """,
        #  ntf   = """
        #            + self.conf_prod["ntf"]
        #            + """,  ntc = """
        #            + self.conf_prod["ntc"]
        #            + """,
        #  nstlim= """
        #            + nstlim
        #            + """, dt= """
        #            + self.conf_prod["dt"]
        #            + """,
        #  cut   = """
        #            + self.conf_prod["cut"]
        #            + """,
        #  temp0 = """
        #            + self.conf_prod["temp0"]
        #            + """,
        #  ntpr  = """
        #            + ntpr
        #            + """, ntwx = """
        #            + self.conf_prod["ntwx"]
        #            + """,
        #  ntt   = """
        #            + self.conf_prod["ntt"]
        #            + """, gamma_ln = """
        #            + self.conf_prod["gamma_ln"]
        #            + """,
        #  ntb   = 2,  ntp = 1,
        #  iwrap = """
        #            + self.conf_prod["iwarp"]
        #            + """,
        #  ig    = -1,
        #"""+ntr_line+nmropt_line+''' /
        #'''+DISANG_tail
        #write
        with open(o_path, 'w') as of:
            of.write(conf_str)
        return o_path


    def _build_MD_rs(self, step, o_path):
        '''
        Generate a file for DISANG restraint. Get parameters from self.conf_step.
        '''
        rs_str = ''
        rs_data_step = self.__dict__['conf_' + step]['rs_constraints']
        for rest_data in rs_data_step:
            rs_str = rs_str + '''  &rst
   iat=  ''' + ','.join(
                rest_data['iat']
            ) + ''', r1= ''' + rest_data['r1'] + ''', r2= ''' + rest_data[
                'r2'] + ''', r3= ''' + rest_data[
                    'r3'] + ''', r4= ''' + rest_data['r4'] + ''',
   rk2=''' + rest_data['rk2'] + ''', rk3=''' + rest_data[
                        'rk3'] + ''', ir6=''' + rest_data[
                            'ir6'] + ''', ialtd=''' + rest_data['ialtd'] + ''',
  &end
'''
        #write
        with open(o_path, 'w') as of:
            of.write(rs_str)
        return o_path

    def reset_MD_conf(self):
        """
        reset MD configuration of current object to default. 
        """
        self.conf_min = Config.Amber.conf_min
        self.conf_heat = Config.Amber.conf_heat
        self.conf_equi = Config.Amber.conf_equi
        self.conf_prod = Config.Amber.conf_prod

    def show_MD_conf(self):
        """
        Show MD configuration of current object. 
        """
        print("Min :     " + repr(self.conf_min))
        print("Heat:     " + repr(self.conf_heat))
        print("Equi:     " + repr(self.conf_equi))
        print("Prod:     " + repr(self.conf_prod))


    def nc2mdcrd(self,
                 o_path="",
                 point=None,
                 start=1,
                 end=-1,
                 step=1,
                 engine="cpptraj"):
        """
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
        if self.nc == None:
            raise Exception(
                "No nc file found. Please assign self.nc or run PDBMD first")
        else:
            if o_path == "":
                o_path = self.nc[:-2] + "mdcrd"
            if end == -1:
                end = "last"
            if point != None:
                all_p = int(self.conf_prod["nstlim"]) / int(
                    self.conf_prod["ntwx"])
                step = int(all_p / point)

            if engine not in ["pytraj", "cpptraj"]:
                raise Exception("engine: pytraj or cpptraj")

            if engine == "pytraj":
                pass

            if engine == "cpptraj":
                cpp_in_path = self.cache_path + "/cpptraj_nc2mdcrd.in"
                cpp_out_path = self.cache_path + "/cpptraj_nc2mdcrd.out"
                with open(cpp_in_path, "w") as of:
                    of.write("parm " + self.prmtop_path + line_feed)
                    of.write("trajin " + self.nc + " " + str(start) + " " +
                             end + " " + str(step) + line_feed)
                    of.write("trajout " + o_path + line_feed)
                    of.write("run" + line_feed)
                    of.write("quit" + line_feed)
                os.system("cpptraj -i " + cpp_in_path + " > " + cpp_out_path)

        self.mdcrd = o_path
        return o_path


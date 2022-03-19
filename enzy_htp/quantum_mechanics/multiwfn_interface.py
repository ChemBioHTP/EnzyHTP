#TODO
class MultiwfnInterface:
    @classmethod
    def get_bond_dipole(cls, qm_fch_paths, a1, a2, prog="Multiwfn"):
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
        Dipoles = []

        if prog == "Multiwfn":
            # self.init_Multiwfn()
            ref_path = qm_fch_paths[0]
            mltwfn_in_path = (ref_path[:-(len(ref_path.split(".")[-1]) + 1)] +
                              "_dipole.in")

            with open(mltwfn_in_path, "w") as of:
                of.write("19" + line_feed)
                of.write("-8" + line_feed)
                of.write("1" + line_feed)
                of.write("y" + line_feed)
                of.write("q" + line_feed)

            bond_id_pattern = r"\( *([0-9]+)[A-Z][A-z]? *- *([0-9]+)[A-Z][A-z]? *\)"
            bond_data_pattern = (
                r"X\/Y\/Z: *([0-9\.\-]+) *([0-9\.\-]+) *([0-9\.\-]+) *Norm: *([0-9\.]+)"
            )

            for fchk in qm_fch_paths:
                # get a1->a2 vector from .out (update to using fchk TODO)
                G_out_path = fchk[:-len(fchk.split(".")[-1])] + "out"
                with open(G_out_path) as f0:
                    coord_flag = 0
                    skip_flag = 0
                    for line0 in f0:
                        if 'Input orientation' in line0:
                            coord_flag = 1
                            continue
                        if coord_flag:
                            # skip first 4 lines
                            if skip_flag <= 3:
                                skip_flag += 1
                                continue
                            if "-------------------------" in line0:
                                break
                            l_p = line0.strip().split()
                            if str(a1) == l_p[0]:
                                coord_a1 = np.array(
                                    (float(l_p[3]), float(l_p[4]),
                                     float(l_p[5])))
                            if str(a2) == l_p[0]:
                                coord_a2 = np.array(
                                    (float(l_p[3]), float(l_p[4]),
                                     float(l_p[5])))
                Bond_vec = coord_a2 - coord_a1

                # Run Multiwfn
                mltwfn_out_path = fchk[:-len(fchk.split(".")[-1])] + "dip"
                if Config.debug >= 2:
                    print("Running: " + Config.Multiwfn.exe + " " + fchk +
                          " < " + mltwfn_in_path)
                run(
                    Config.Multiwfn.exe + " " + fchk + " < " + mltwfn_in_path,
                    check=True,
                    text=True,
                    shell=True,
                    capture_output=True,
                )
                run(
                    "mv LMOdip.txt " + mltwfn_out_path,
                    check=True,
                    text=True,
                    shell=True,
                    capture_output=True,
                )
                run(
                    "rm LMOcen.txt new.fch",
                    check=True,
                    text=True,
                    shell=True,
                    capture_output=True,
                )

                # get dipole
                with open(mltwfn_out_path) as f:
                    read_flag = 0
                    for line in f:
                        if line.strip(
                        ) == "Two-center bond dipole moments (a.u.):":
                            read_flag = 1
                            continue
                        if read_flag:
                            if "Sum" in line:
                                raise Exception("Cannot find bond:" + str(a1) +
                                                "-" + str(a2) + line_feed)
                            Bond_id = re.search(bond_id_pattern, line).groups()
                            # find target bond
                            if str(a1) in Bond_id and str(a2) in Bond_id:
                                Bond_data = re.search(bond_data_pattern,
                                                      line).groups()
                                dipole_vec = (
                                    float(Bond_data[0]),
                                    float(Bond_data[1]),
                                    float(Bond_data[2]),
                                )
                                # determine sign
                                if np.dot(np.array(dipole_vec), Bond_vec) > 0:
                                    dipole_norm_signed = float(Bond_data[3])
                                else:
                                    dipole_norm_signed = -float(Bond_data[3])
                                break
                Dipoles.append((dipole_norm_signed, dipole_vec))

        return Dipoles

    @classmethod
    def init_Multiwfn(cls, n_cores=None):
        """
        initiate Multiwfn with settings in Config
        """
        # set nthreads
        if n_cores == None:
            n_cores = str(Config.n_cores)
        if Config.debug >= 1:
            print("Running: " + "sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  " +
                  n_cores + "/' " + Config.Multiwfn.DIR + "/settings.ini")
        run(
            "sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  " + n_cores + "/' " +
            Config.Multiwfn.DIR + "/settings.ini",
            check=True,
            text=True,
            shell=True,
            capture_output=True,
        )


def get_fchk(self, keep_chk=0):
    """
    transfer Gaussian chk files to fchk files using formchk
    ----------
    keep_chk: if not delete original chk file (default: 0)
    """
    # san check
    if len(self.qm_cluster_chk) == 0:
        raise Exception("No chk file in self.qm_cluster_chk.")

    # formchk
    fchk_paths = []
    for chk in self.qm_cluster_chk:
        fchk = chk[:-3] + "fchk"
        if Config.debug > 1:
            print("running: " + "formchk " + chk + " " + fchk)
        run(
            "formchk " + chk + " " + fchk,
            check=True,
            text=True,
            shell=True,
            capture_output=True,
        )
        fchk_paths.append(fchk)
        # keep chk
        if not keep_chk:
            if Config.debug > 1:
                print("removing: " + chk)
            os.remove(chk)

    self.qm_cluster_fchk = fchk_paths
    return self.qm_cluster_fchk

"""calculate DSI of given traj file
TODO make this into an API under analysis
NOTE(2025.8): An API of DSI is made in PR #234"""
from typing import List

from enzy_htp import interface
import enzy_htp.core.file_system as fs
import numpy as np

def dsi(traj_file: str, top_file: str,
        domain1: tuple, domain2: tuple, 
        int_file: str = None):
    """calculate a list of DSI for each frame from the traj_file"""
    domain1_pattern = f":{domain1[0]}-{domain1[1]}&!@H="
    domain2_pattern = f":{domain2[0]}-{domain2[1]}&!@H="
    if not int_file:
        int_file = fs.get_valid_temp_name("dsi.dat")
    contents: List[str] = [
        f"parm {top_file}",
        f"trajin {traj_file}",
        f"distance d_CD-CBM {domain1_pattern} {domain2_pattern} out {int_file} geom",
        f"radgyr Rg_CD {domain1_pattern} out {int_file} nomax",
        f"radgyr Rg_CBM {domain2_pattern} out {int_file} nomax",
        "run",
        "quit",
    ]
    contents = "\n".join(contents)
    interface.amber.run_cpptraj(contents)

    result = []
    with open(int_file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            index, distance, rd1, rd2 = line.strip().split()
            result.append(float(distance) - float(rd1) - float(rd2))

    fs.clean_temp_file_n_dir(int_file)

    return result

def main():
    linker_list = [
        "GSGDGGGNDGGEGGL",
        "DGPAPEPVKHPIDHVG",
        "FKAQPDLAEAAATTTENP",
        "GTLSP",
        "ENTRGFHG",
        "PNIATG",
        "ESPVDSEQRGENDL",
        "LSAYERSCGIPDRM",
        "NETTNK",
        "FWGMASSSY",
        "GSAGSAAGSGEF",
        "MEQAM",
        "VKGKEGDEEEE",
        "LDSGA",
        "SDP",
        "DYGNSPLHRFKKPGSKNFQNIFPPSAT",
        "ETGLNAYLPGLAGKE",
        "KTWENVNAQ",
        "YDRPDR",
        "TNGRL",
        "NLEEHLGKLN",
        "RAETIDDIR",
        "AKLKQKTEQLQDRIAG",
        "LEEVPSVGVNKNIFL",
        "ITAKDE",
        "LNRLDRL",
    ]
    domain_list = [[(1, 418), (419+len(i), 515+len(i))] for i in linker_list]

    dsi_values = dsi("samples.xtc", "topology.pdb", *domain_list[-2])
    fs.write_lines("dsi.csv", map(str, dsi_values))
    print(np.mean(dsi_values))

if __name__ == "__main__":
    main()

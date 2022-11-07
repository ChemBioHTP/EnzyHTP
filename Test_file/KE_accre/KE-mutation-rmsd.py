import pathlib
from glob import glob
from subprocess import run
import os

from Class_PDB import PDB

mask = ":9,11,48,50,101,128,201,202,222&!@H=" # 5A
mask_sasa = ":9,11,48,50,101,128,201,202,222"
mask_pro = ":1-253"
mask_sub = ":254"

def main():
    for i in range(6):
        target_groups = glob(f"/data/yang_lab/shaoqz/KE-DE/R5/group_{i}*")
        for mut_group in target_groups:

            with open(f"{mut_group}/Mutation.dat") as f:
                lines = f.readlines()
                insert_mapper = {}
                for j, line in enumerate(lines):
                    if "traj" in line:
                        traj_path = pathlib.Path(f"{mut_group}/{eval(lines[j+1].strip())}")
                        insert_idx = j+2
                        prmtop_path = list(traj_path.parent.parent.glob("*prmtop"))[0]
                        
                        rmsd_value = PDB.get_rmsd(str(prmtop_path), str(traj_path), mask)
                        sasa_value = PDB.get_sasa_ratio(str(prmtop_path), str(traj_path), 
                                                        mask_pro, mask_sasa, mask_sub)
                        insert_mapper[insert_idx] = f"---rmsd_avg---\n{repr(rmsd_value)}\n---sasa_avg---\n{sasa_value}\n"
            
            with open(f"./Mutation_{i}.dat", "a") as of:
                add_line = 0
                for li, value in insert_mapper.items():
                    lines.insert(li+add_line, value)
                    add_line += len(value.strip().split("\n"))
                of.writelines(lines)
            break
        break


if __name__ == "__main__":
    main()


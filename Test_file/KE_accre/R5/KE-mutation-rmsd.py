import pathlib
from glob import glob
import os

from Class_PDB import PDB

mask = ":9,11,48,50,101,128,201,202,222&!@H=" # 5A
mask_sasa = ":9,11,48,50,101,128,201,202,222"
mask_pro = ":1-253"
mask_sub = ":254"

def get_rmsd_sasa_all():
    for i in range(6):
        target_groups = glob(f"/data/yang_lab/shaoqz/KE-DE/R5/group_{i}*")
        if i == 0:
            target_groups += glob(f"/data/yang_lab/shaoqz/KE-DE/R5/group_sele*")
        for mut_group in target_groups:
            if not os.path.exists(f"{mut_group}/Mutation.dat"):
                print(f"no data file under {mut_group}")
                continue
            with open(f"{mut_group}/Mutation.dat") as f:
                lines = f.readlines()
                insert_mapper = {}
                for j, line in enumerate(lines):
                    if "TAG" in line:
                        muta_flag_str = "_".join(["".join(x) for x in eval(lines[j+1].strip())])
                    if "traj" in line:
                        traj_path = pathlib.Path(f"{mut_group}/{eval(lines[j+1].strip())}")
                        insert_idx = j+2
                        prmtop_path = list(traj_path.parent.parent.glob(f"*{muta_flag_str}*prmtop"))[0]
                        
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

def get_rmsd_sasa_part():
    mut_group = "/data/yang_lab/shaoqz/KE-DE/R5/group_2_1"
    if not os.path.exists(f"{mut_group}/Mutation.dat"):
        Exception(f"no data file under {mut_group}")
    with open(f"{mut_group}/Mutation.dat") as f:
        lines = f.readlines()
        insert_mapper = {}
        for j, line in enumerate(lines):
            if "TAG" in line:
                muta_flag_str = "_".join(["".join(x) for x in eval(lines[j+1].strip())])
            if "traj" in line:
                traj_path = pathlib.Path(f"{mut_group}/{eval(lines[j+1].strip())}")
                insert_idx = j+2
                prmtop_path = list(traj_path.parent.parent.glob(f"*{muta_flag_str}*prmtop"))[0]
                
                rmsd_value = PDB.get_rmsd(str(prmtop_path), str(traj_path), mask)
                sasa_value = PDB.get_sasa_ratio(str(prmtop_path), str(traj_path), 
                                                mask_pro, mask_sasa, mask_sub)
                insert_mapper[insert_idx] = f"---rmsd_avg---\n{repr(rmsd_value)}\n---sasa_avg---\n{sasa_value}\n"
    
    with open(f"./Mutation_add.dat", "a") as of:
        add_line = 0
        for li, value in insert_mapper.items():
            lines.insert(li+add_line, value)
            add_line += len(value.strip().split("\n"))
        of.writelines(lines)

def replace_EF_dipole_data():
    for i in range(6):
        target_groups = glob(f"/data/yang_lab/shaoqz/KE-DE/R5/group_{i}*")
        if i == 0:
            target_groups += glob(f"/data/yang_lab/shaoqz/KE-DE/R5/group_sele*")
        new_data_mapper = {}
        # get new data mapping with muta_flag_line
        for mut_group in target_groups:
            data_file_path = f"{mut_group}/Mutation.dat"
            if not os.path.exists(data_file_path):
                print(f"no data file under {mut_group}")
                continue
            with open(data_file_path) as f:
                lines = f.readlines()
                insert_mapper = {}
                for j, line in enumerate(lines):
                    if "TAG" in line:
                        muta_flag_line = lines[j+1]
                    if "traj" in line:
                        new_data_file_path = list(pathlib.Path(f"{mut_group}/{eval(lines[j+1].strip())}").parent.parent.glob("EF_QM.dat"))[0]
                        with open(new_data_file_path) as f1:
                            new_EF_dipole_lines = f1.readlines()[2:]
                        new_data_mapper[muta_flag_line] = new_EF_dipole_lines
        with open(f"./Mutation_{i}.dat") as f2:
            f2_lines = f2.readlines()
            of_lines = []
            k = 0
            while k < len(f2_lines):
                if "TAG" in f2_lines[k]:
                    new_data = new_data_mapper[f2_lines[k+1]]
                    of_lines.append(f2_lines[k])
                    of_lines.append(f2_lines[k+1])
                    of_lines.extend(new_data)
                    k += len(new_data) +2
                    continue
                of_lines.append(f2_lines[k])
                k += 1
            with open(f"./new/Mutation_{i}.dat", "w") as of:
                of.writelines(of_lines)
                                    
def main():
    # get_rmsd_sasa_all()
    # get_rmsd_sasa_part()
    replace_EF_dipole_data()

if __name__ == "__main__":
    main()

    
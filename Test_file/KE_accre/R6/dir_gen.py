from glob import glob
import pathlib
from subprocess import run
import shutil

from helper import mkdir

def generate_input_dirs():
    start_stru = glob("start_structure/*pdb")

    for i, stru_path in enumerate(start_stru):
        for j in range(5):
            dir_name = f"./group_{i}_{j}"
            mkdir(dir_name)
            # move source pdb
            pdb_name = stru_path.removeprefix("start_structure/")
            shutil.copy(stru_path, dir_name)
            # move script
            script_path = shutil.copy("./KE-metrics_complete.py", dir_name)
            run(f'sed -i "s/XXX/{pdb_name}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
            run(f'sed -i "s/YYY/{2}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
            # move submission script
            shutil.copy("./sub_enzy_htp.cmd", dir_name)
            mkdir(f"{dir_name}/ligands")
            shutil.copy("./ligands/ligand_H5J.frcmod", f"{dir_name}/ligands/")
            shutil.copy("./ligands/ligand_H5J.prepin", f"{dir_name}/ligands/")
    with open("group_mapper.mapper", "w") as of:
        of.write(repr(start_stru))

def read_group_mapper():
    with open("group_mapper.mapper") as f:
        group_stru_list = eval(f.read())
    for stru_path in group_stru_list:
        muta_flags = pathlib.Path(stru_path).stem.removeprefix("KE-07_").split("_")
        muta_flags = [flag[:1]+flag[2:] for flag in muta_flags]
        print(muta_flags)
def main():
    read_group_mapper()

if __name__ == "__main__":
    main()



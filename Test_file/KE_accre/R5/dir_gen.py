from glob import glob
from subprocess import run
import shutil

from helper import mkdir

start_stru = glob("start_structure/*pdb")

for i, stru_path in enumerate(start_stru):
    for j in range(5):
        dir_name = f"./group_{i}_{j}"
        mkdir(dir_name)
        # move source pdb
        pdb_name = stru_path.removeprefix("start_structure/")
        shutil.copy(stru_path, dir_name)
        # move script
        script_path = shutil.copy("./KE-metrics.py", dir_name)
        run(f'sed -i "s/XXX/{pdb_name}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
        run(f'sed -i "s/YYY/{2}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
        # move submission script
        shutil.copy("./sub_enzy_htp.cmd", dir_name)




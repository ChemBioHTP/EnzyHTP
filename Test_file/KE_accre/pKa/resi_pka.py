import os
from subprocess import run
from Class_PDB import PDB
from glob import glob

data_dir = "/scratch/jiany37/KE07-R7-2_CHdirec_mut/5_MD/"
target_mutants = {
    "E24K-S29K-K162F-R163F" : [-3.099821946],
    "E24G-S29R-K162A-R163Y" : [-2.992261664],
    "E24R-S29R-K162Q-R163F" : [-2.940158084],
    "E24N-S29H-K162Y-R163Y" : [-2.919226723],
    "E24Q-S29K-K162Y-R163Y" : [-2.919066654],
    "E24Q-S29H-K162Q-R163F" : [-2.744995623],
    "E24Q-S29K-K162Q-R163F" : [-2.242331157],
    "E24F-S29R-K162W-R163L" : [-3.302921812],
    "E24L-S29R-K162C-R163L" : [-3.124140161],
    "E24I-S29R-K162C-R163L" : [-3.063535635],
    "E24M-S29H-K162L-R163L" : [-2.952885940],
    "E24M-S29H-K162W-R163L" : [-2.734547929],
    "E24M-S29H-K162C-R163L" : [-2.684644251],
    "E24L-S29R-K162V-R163L" : [-2.679997883],
    "E24M-S29H-K162L-R163Y" : [-2.626222923],
    "WT": [0.0]
}
temp_pdb_path = "./temp_prod.pdb"


def main():
    for mutant in target_mutants:
        if mutant == "WT":
            rst_path = "/scratch/stullsl/KE07Simulation/WT/WT1/MD/prod.rst"
            prmtop_path = glob("/scratch/stullsl/KE07Simulation/WT/WT1/*prmtop")[0]
        else:
            rst_path = f"{data_dir}{mutant}/MD/prod.rst"
            prmtop_path = glob(f"{data_dir}{mutant}/*prmtop")[0]
        run(f"ambpdb -p {prmtop_path} -c {rst_path} > {temp_pdb_path}",
            check=True, text=True, shell=True, capture_output=True)
        # restore HIS and remove substrate, water and Na+
        run(rf"sed -i 's/HIE/HIS/g;s/HID/HIS/g;s/HIP/HIS/g;s/.*H5J.*//g;s/.*WAT.*//g;s/.*Na+.*//g;/^[[:space:]]*$/d' {temp_pdb_path}",
            check=True, text=True, shell=True, capture_output=True)
        pdb_obj = PDB(temp_pdb_path)
        target_pka = pdb_obj.get_residue_pka(":101")[0]
        target_mutants[mutant].append(target_pka)
        # clean up
        os.remove(temp_pdb_path)

    with open("pka_E101_out.csv", "w") as of:
        for mutant, m_data in target_mutants.items():
            of.write(f"{mutant},{','.join(str(x) for x in m_data)}{os.linesep}")

if __name__ == "__main__":
    main()


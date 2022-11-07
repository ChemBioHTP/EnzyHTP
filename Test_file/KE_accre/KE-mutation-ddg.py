import copy
from glob import glob
import re

from core.clusters.accre import Accre
from Class_PDB import PDB
from Class_Conf import Config
from helper import write_data

data_output_path = './Mutation.dat'

# index
variant_idx = 0

def get_unique_mutant(mutaflags: list) -> list:
    unique_mutation = {}
    for mutaflag in mutaflags:
        position = mutaflag[1]
        if position in unique_mutation.keys():
            print(f"found duplicated mutation in {position}, using {mutaflag}")
        unique_mutation[position] = mutaflag
    return unique_mutation.values()
        

def main():
    start_pdb = PDB('./KE-07_aH_IA7D_KA146E_GA202R_NA224D.pdb')
    # calculation mutants
    initial_rel_mutant_list = [
        "D7Q I199V R202G F227L", "D7Q F229S", "",
        "D7T E146T R202G D224N I199Q F86L I173V L176D F227L", "D7V R202G I199Q K162P I173A L176I F229S",
        "D7S E146K K19E"]
    muta_groups = []
    for i, init_mut in enumerate(initial_rel_mutant_list):
        init_mutaflags = []
        if init_mut:
            for mut in init_mut.split(" "):
                r_mut = (mut[0], mut[1:-1], mut[-1])
                init_mutaflags.append(r_mut)
        with open(f"Mutation_{i}.dat") as f:
            add_muts = re.findall(r"TAG===\n(.+?)\n", f.read())
            for add_mut in add_muts:
                overall_mut = copy.deepcopy(init_mutaflags)
                add_mut = eval(add_mut)
                for add_mut_i in add_mut:
                    add_mut_i = add_mut_i[0:1]+add_mut_i[2:]
                    overall_mut = overall_mut + [add_mut_i]
                unique_mutaflags = get_unique_mutant(overall_mut)
                muta_groups.append(tuple(map(lambda x: " ".join(x), unique_mutaflags)))

    # relax
    relaxed_pdb = start_pdb.relax_with_rosetta(
        rosetta_home='/data/yang_lab/shaoqz/software/Rosetta313/main/',
        cluster=Accre(),
        period=30)
    # run ddg
    ddg_results = start_pdb.get_rosetta_ddg(
        rosetta_home='/data/yang_lab/Common_Software/Rosetta3.9/main/',
        muta_groups=muta_groups,
        relaxed_pdb=relaxed_pdb,
        job_array_size=100,
        cluster=Accre(),
        period=30)
    with open("ddg_results.dat", "w") as of:
        of.write(repr(ddg_results))

if __name__ == "__main__":
    main()

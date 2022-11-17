import copy
from glob import glob
import re
import sys


from core.clusters.accre import Accre
from Class_PDB import PDB

data_output_path = "ddg_results.dat"
initial_mutant_list = [
    ['I7D', 'K146E', 'G202R', 'N224D', 'I151T'],
    ['I7Q', 'K146E', 'G202R', 'N224D', 'F229S', 'N103Q'],
    ['I7S', 'K19E', 'G202R', 'N224D', 'V190M'],
    ['I7S', 'K19E', 'G202R', 'N224D', 'E239M'],
    ['I7T', 'K146T', 'I199Q', 'F86L', 'I173V', 'L176D', 'F227L', 'E64W'],
    ['I7T', 'K146T', 'I199Q', 'F86L', 'I173V', 'L176D', 'F227L', 'E231S'],
    ['I7D', 'G202R', 'N224D', 'V12M'],
    ['I7Q', 'K146E', 'G202R', 'N224D', 'F229S', 'K132S'],
    ['I7S', 'K19E', 'G202R', 'N224D', 'S144H'],
    ['I7D', 'K146E', 'G202R', 'N224D', 'I102L'],
    ['I7Q', 'K146E', 'G202R', 'N224D', 'F229S', 'P197W'],
    ['I7D', 'K146E', 'G202R', 'N224D', 'L47I']]

def decode_mutaflag(mutaflag: str) -> tuple:
    old_resi = mutaflag[0]
    new_resi = mutaflag[-1]
    position = mutaflag[1:-1]
    return old_resi, new_resi ,position

def get_unique_mutant(mutaflags: list) -> list:
    unique_mutation = {}
    for mutaflag in mutaflags:
        position = mutaflag[1]
        if position in unique_mutation:
            old_flag = unique_mutation[position]
            if mutaflag[0] == old_flag[-1]:
                if old_flag[0] != mutaflag[-1]:
                    print(f"found additional mutation in {position} of {mutaflags}, update it")
                    unique_mutation[position] = (old_flag[0], position, mutaflag[-1])
                else:
                    print(f"mutation back in {position} of {mutaflags}, delete both")
                    del unique_mutation[position]
            else:
                raise Exception(f"inconsistant relative residue in position {position}: {mutaflags}")
            continue
        unique_mutation[position] = mutaflag
    return unique_mutation.values()

def calculation_relative_mutant(mut_flags: list, ref_mut_flags: list) -> list:
    """
    ['A' '11' 'B', 'A' '12' 'C'] ref: ['A' '11' 'D', 'A' '13' 'D'] -> ['D' '11' 'B', 'A' '12' 'B', 'D' '13' 'A']
    """
    ref_mut_mapper = {}
    result = []
    for r_flag in ref_mut_flags:
        ref_mut_mapper[r_flag[1]] = r_flag
    for flag in mut_flags:
        position = flag[1]
        if position in ref_mut_mapper:
            ref_flag = ref_mut_mapper[position]
            if flag[0] == ref_flag[0]:
                if flag[-1] == ref_flag[-1]:
                    #print(f"found same mutation {flag}, detele it in relative representation")
                    del ref_mut_mapper[position]
                    continue
                #print(f"found additional mutation {flag}, update it in relative representation")
                flag = (ref_flag[-1], position, flag[-1])
            else:
                raise Exception(f"inconsistant WT residue in position {position}: {mut_flags}")
            del ref_mut_mapper[position]
        result.append(flag)
    
    for position_left, ref_flag in ref_mut_mapper.items():
        new_flag = (ref_flag[-1], ref_flag[1], ref_flag[0])
        result.append(new_flag)
    return result

def get_ddg_part():
    start_pdb = PDB('./KE-07_aH_IA7D_KA146E_GA202R_NA224D.pdb')
    # calculation mutants
    init_mut = "D7S E146K K19E"
    muta_groups = []
    init_mutaflags = []
    if init_mut:
        for mut in init_mut.split(" "):
            r_mut = (mut[0], mut[1:-1], mut[-1])
            init_mutaflags.append(r_mut)
    with open(f"Mutation_add.dat") as f:
        add_muts = re.findall(r"TAG===\n(.+?)\n", f.read())
        for add_mut in add_muts:
            overall_mut = copy.deepcopy(init_mutaflags)
            add_mut = eval(add_mut)
            for add_mut_i in add_mut:
                add_mut_i = add_mut_i[0:1]+add_mut_i[2:]
                overall_mut = overall_mut + [add_mut_i]
            unique_mutaflags = get_unique_mutant(overall_mut)
            muta_groups.append(tuple(map(lambda x: " ".join(x), unique_mutaflags)))

    # run ddg
    ddg_results = start_pdb.get_rosetta_ddg(
        rosetta_home='/data/yang_lab/Common_Software/Rosetta3.9/main/',
        muta_groups=muta_groups,
        relaxed_pdb='/home/shaoq1/KE-DE/R5/ddg_calc/KE-07_aH_IA7D_KA146E_GA202R_NA224D_relaxed.pdb',
        job_array_size=100,
        cluster=Accre(),
        period=30)
    with open("ddg_result_add.dat", "w") as of:
        of.write(repr(ddg_results))

def get_ddg_all_w_relax():
    start_pdb = PDB('./KE-07_aH_IA7D_KA146E_GA202R_NA224D.pdb')
    start_mutation = ['I7D','K146E','G202R','N224D']
    # calculate relative mutations
    muta_groups = []
    for i, init_mut in enumerate(initial_mutant_list):
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
                rel_mutaflags = calculation_relative_mutant(unique_mutaflags)
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
    with open(data_output_path, "w") as of:
        of.write(repr(ddg_results))

def get_ddg_all_wo_relax():
    start_pdb = PDB('./KE-07_aH_IA7D_KA146E_GA202R_NA224D.pdb')
    relaxed_pdb = './KE-07_aH_IA7D_KA146E_GA202R_NA224D_relaxed.pdb'
    start_mutation = ['I7D','K146E','G202R','N224D']
    start_r_mutaflags = []
    for mut in start_mutation:
        start_r_mutaflags.append((mut[0], mut[1:-1], mut[-1]))

    # calculate relative mutations
    muta_groups = []
    for i, init_mut in enumerate(initial_mutant_list):
        init_mutaflags = []
        if init_mut:
            for mut in init_mut:
                r_mut = (mut[0], mut[1:-1], mut[-1])
                init_mutaflags.append(r_mut)
        with open(f"Mutation_{i}.dat") as f:
            variants = re.findall(r"TAG===\n(.+?)\n", f.read())
            for add_mut in variants:
                unique_add_mut = {}
                overall_mut = copy.deepcopy(init_mutaflags)
                add_mut = eval(add_mut)
                for add_mut_i in add_mut:
                    add_mut_i = add_mut_i[0:1]+add_mut_i[2:]
                    if add_mut_i[1] in unique_add_mut:
                        print(f"found duplicated mutations from EnzyHTP: {add_mut}, later one is used")
                    unique_add_mut[add_mut_i[1]] = add_mut_i
                overall_mut.extend(unique_add_mut.values())
                unique_mutaflags = get_unique_mutant(overall_mut)
                rel_mutaflags = calculation_relative_mutant(unique_mutaflags, start_r_mutaflags)
                muta_groups.append(tuple(map(lambda x: " ".join(x), rel_mutaflags)))
    # run ddg
    ddg_results = start_pdb.get_rosetta_ddg(
        rosetta_home='/data/yang_lab/Common_Software/Rosetta3.9/main/',
        muta_groups=muta_groups,
        relaxed_pdb=relaxed_pdb,
        job_array_size=100,
        cluster=Accre(),
        period=30)
    with open(data_output_path, "w") as of:
        of.write(repr(ddg_results))

def main():
    get_ddg_all_wo_relax()

if __name__ == "__main__":
    main()

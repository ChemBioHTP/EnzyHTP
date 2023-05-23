import copy
import re
import pickle


from core.clusters.accre import Accre
from Class_PDB import PDB

data_output_path = "ddg_results.dat"
resi_idx_mapper = {
    "TynA":{
        "A": -5,
        "B": 715
        },
    "PuO":{
        "A": -1,
        "B": 449 
    }
}

def decode_mutaflag(mutaflag: str) -> tuple:
    old_resi = mutaflag[0]
    new_resi = mutaflag[-1]
    chain_id = mutaflag[1]
    position = mutaflag[2:-1]
    return old_resi, new_resi ,chain_id, position

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

def parse_mutation_to_rosetta(mutations, mutation_mapper):
    result = []
    for mut in mutations:
        old_resi, new_resi, chain_id, position = decode_mutaflag(mut)
        position = str(int(position) + mutation_mapper[chain_id])
        result.append(" ".join((old_resi, position, new_resi)))
    return tuple(result)

def get_ddg_all_w_relax_acpu(pdb_path: str, mutants: list, mutation_mapper: dict):
    start_pdb = PDB(pdb_path)
    relaxed_pdb = f'{pdb_path.removesuffix(".pdb")}_relaxed.pdb'

    muta_groups = []
    for mut in mutants:
        if mut is []:
            continue
        r_mutaflags = parse_mutation_to_rosetta(mut, mutation_mapper)
        muta_groups.append(r_mutaflags)
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

def get_ddg_all_wo_relax_acpu(pdb_path: str, mutants: list, mutation_mapper: dict):
    start_pdb = PDB(pdb_path)
    relaxed_pdb = f'{pdb_path.removesuffix(".pdb")}_relaxed.pdb'

    muta_groups = []
    for mut in mutants:
        if mut is []:
            continue
        r_mutaflags = parse_mutation_to_rosetta(mut, mutation_mapper)
        muta_groups.append(r_mutaflags)
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
    with open("./mutants.pickle", "rb") as f:
        mutants = pickle.load(f)
    get_ddg_all_wo_relax_acpu("./PuOrh_amber_aH.pdb", mutants , resi_idx_mapper["PuO"])


if __name__ == "__main__":
    main()

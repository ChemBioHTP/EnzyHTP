from typing import List
from itertools import product
import pickle

WORK_DIR = "/data/yang_lab/shaoqz/AcPutr_detect/TOT/1st/"

mutation_section_tyna = {
    "T223":"RHKYC0",
    "L225":"RHKYC0",
    "Y381":"RHKC0"
    }

mutation_section_puo = {
    "E324":"RHKYC0",
    "G172":"RHKDESTNQCPAVILMFYW0"
    }

mutation_section_puo_2nd = {
    "E324":"N",
    "G172":"RHKDESTNQCPAVILMFYW0"
    }

def timer(fn):
    from time import perf_counter
    
    def inner(*args, **kwargs):
        start_time = perf_counter()
        to_execute = fn(*args, **kwargs)
        end_time = perf_counter()
        execution_time = end_time - start_time
        print('{0} took {1:.8f}s to execute'.format(fn.__name__, execution_time))
        return to_execute
    
    return inner

def combine_lists(target_list):
    curr_list = next(target_list, None)
    if not curr_list:
        return [[]]
    next_list = combine_lists(target_list)
    return [[x]+y for x in curr_list for y in next_list]

def combine_lists_wo_empty(target_list):
    curr_list = next(target_list, None)
    if not curr_list:
        return [[]]
    next_list = combine_lists_wo_empty(target_list)
    return [[x]+y if x != "" else y for x in curr_list for y in next_list]

@timer
def combine_mutations(mutation_mapper: dict) -> List[str]:
    """generate combinatorial mutation space"""
    result = []
    for posi, targets in mutation_mapper.items():
        posi_muts = []
        holder = []
        for target in targets:
            mut = posi+target
            posi_muts.append([mut])
 
        for old_mut in result:
            for mut in posi_muts:
                holder.append(old_mut + mut)
        if result:
            result = holder
        else:
            result = posi_muts
    return result

@timer
def combine_mutations_2(mutation_mapper: dict) -> List[str]:
    """generate combinatorial mutation space"""
    result = []
    for posi, targets in mutation_mapper.items():
        posi_muts = []
        for target in targets:
            mut = posi+target
            posi_muts.append(mut)
        result.append(posi_muts)
    result = combine_lists(iter(result))
    return result

@timer
def combine_mutations_3(mutation_mapper: dict) -> List[str]:
    """generate combinatorial mutation space"""
    result = []
    for posi, targets in mutation_mapper.items():
        posi_muts = []
        for target in targets:
            mut = posi+target
            posi_muts.append(mut)
        posi_muts.append("")
        result.append(posi_muts)
    result = product(*result)
    return result

@timer
def combine_mutations_4(mutation_mapper: dict) -> List[str]:
    """generate combinatorial mutation space. including not mutate a position."""
    result = []
    for posi, targets in mutation_mapper.items():
        posi_muts = []
        for target in targets:
            if target == "0":
                mut = ""
            else:
                mut = posi+target
            posi_muts.append(mut)
        result.append(posi_muts)
    result = combine_lists_wo_empty(iter(result))
    return result


def sync_over_chain(mutant_list: List[List[str]], chain_sync_list: List[str]) -> List[List[str]]:
    """sync mutation of ploymer over chains indicated in the chain_sync_list"""
    result = []
    for mutant in mutant_list:
        new_mutant = []
        for mutation in mutant:
            for chain_id in chain_sync_list:
                new_mutant.append(mutation[0]+chain_id+mutation[1:])
        result.append(new_mutant)
    return result

def mutant_gen(mutation_mapper, chain_sync_list):
    return sync_over_chain(combine_mutations_4(mutation_mapper), chain_sync_list)

def main():
    # print(*mutant_gen(mutation_section_tyna, ["A","B"]), sep='\n')
    # with open(f"{WORK_DIR}TynA/mutants.pickle", "wb") as of:
    #     pickle.dump(mutant_gen(mutation_section_tyna, ["A", "B"]), of, 0)
    # with open(f"{WORK_DIR}PuO/mutants.pickle", "wb") as of:
    #     pickle.dump(mutant_gen(mutation_section_puo, ["A", "B"]), of, 0)   
    with open(f"{WORK_DIR}PuO/mutants_2nd.pickle", "wb") as of:
        pickle.dump(mutant_gen(mutation_section_puo_2nd, ["A", "B"]), of, 0)   
if __name__ == "__main__":
    main()

import os
import pickle
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import repeat
from scipy.stats import probplot

#Set project path here
proj_path='./'
csv_file = f'{proj_path}PuO_binding.csv'
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
substrate_list = ["acp","cad","hex","lyl","put","spe"]

def convert_index(chain_id, index, enzy):
    new_index = index-resi_idx_mapper[enzy][chain_id]
    return new_index

def merge_multi_chain_mutation(mutant):
    uni_index_mutant = {}
    for mut in mutant:
        uni_index_mutant[mut[2:-1]] = mut[0]+mut[2:]
    return list(uni_index_mutant.values())


def data_extraction_multi_sub(substrate_list: list, data_dir: str) -> list:
    """
    take inital mutant list and the data dir, extract data for each inital mutant
    """
    data = {}
    mutant_data_pattern = r"===TAG===\n(.+?)\n((?:---[A-z, ]+---\n.+?\n)+)"
    metric_data_pattern = r"---([A-z, ]+)---\n(.+?)\n"

    for sub_name in substrate_list:
        sub_data = []
        with open(f'{data_dir}/{sub_name}_result.dat') as f:
            mutant_data_list = re.findall(mutant_data_pattern, f.read())
            for mutant_data in mutant_data_list:
                # prepare muta flag
                muta_flag_list = eval(mutant_data[0])
                muta_flag_list = list(map(
                    lambda x: "".join(
                        [x[0],
                         x[1],
                         str(convert_index(x[1],int(x[2]),"PuO")),
                         x[3]]),
                    muta_flag_list))
    
                # prepare metric data
                metric_data_list = re.findall(metric_data_pattern, mutant_data[1])
                metric_data_mapper = {}
                for metric_data in metric_data_list:
                    metric_data_mapper[metric_data[0]] = eval(metric_data[1])
                sub_data.append((muta_flag_list, metric_data_mapper))
        data[sub_name] = sub_data
    return data

data = data_extraction_multi_sub(substrate_list, proj_path)
    
# Output
#=========
#---csv---
with open("./mutants.pickle", "rb") as f:
    mutant_list = pickle.load(f)

with open(csv_file,'w') as of:
    # title
    of.write("mutant,")
    for sub_name, sub_data in data.items():
        m_name_0, m_data_0 = sub_data[0]
        of.write(f"{','.join(map(lambda x: f'{sub_name}_{x}', m_data_0.keys()))},")
    of.write(os.linesep)
    # data_line
    for mutant in mutant_list:
        data_line = " ".join(merge_multi_chain_mutation(mutant))+","
        for sub_name, sub_data in data.items():
            sub_data_str = ",,"
            for m_name, m_data in sub_data:
                if set(m_name) == set(mutant):
                    sub_data_str = ','.join(map(lambda x: str(x), m_data.values()))+","
                    break
            data_line += sub_data_str
        of.write(f"{data_line}{os.linesep}")

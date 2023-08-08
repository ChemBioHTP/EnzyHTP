import csv
import os
import pickle
import re
import sys
import numpy as np
from numpy.core.fromnumeric import repeat
from scipy.stats import probplot

#Set project path here
proj_path='./'
csv_file = f'{proj_path}sasa_conv_test.csv'

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
                muta_flag = "WT"
                # prepare metric data
                metric_data_list = re.findall(metric_data_pattern, mutant_data[1])
                metric_data_mapper = {}
                for metric_data in metric_data_list:
                    metric_data_mapper[metric_data[0]] = eval(metric_data[1])
                sub_data.append(metric_data_mapper)
        data[sub_name] = sub_data
    return data

data = data_extraction_multi_sub(["total"], proj_path)
    
# Output
#=========
#---csv---
w_data = [list(i.values()) for i in data["total"]]
w_data = [(np.mean(i[0]), i[1]) for i in w_data]
with open(csv_file,'w', newline="") as of:
    writer = csv.writer(of)
    writer.writerow(["ef", "SPI"])
    writer.writerows(w_data)

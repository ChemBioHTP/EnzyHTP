import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import repeat
from scipy.stats import probplot

#Set project path here
proj_path='./'
csv_file = f'{proj_path}R5.csv'
initial_mutant_list = [
    "I7Q K146E I199V N224D F227L", "I7Q K146E G202R N224D F229S",
    "I7D K146E G202R N224D", "I7T K146T I199Q F86L I173V L176D F227L",
    "I7V K146E I199Q N224D L162P I173A L176I F229S", "I7S K19E G202R N224D"]


def get_unique_mutant(mutaflag_str: list) -> list:
    """
    'A11B A12B A11X' -> 'A12B A11X'
    """
    unique_mutation = {}
    for mutaflag in mutaflag_str.split(" "):
        position = int(mutaflag[1:-1])
        if position in unique_mutation.keys():
            print(f"found duplicated mutation in {position}, using {mutaflag}")
        unique_mutation[position] = mutaflag
    return " ".join(list(unique_mutation.values()))

def data_extraction_multiple_initial(inital_mutant_list: list, data_dir: str) -> list:
    """
    take inital mutant list and the data dir, extract data for each inital mutant
    """
    data = []
    mutant_data_pattern = r"===TAG===\n(.+?)\n((?:---[A-z, ]+---\n.+?\n)+)"
    metric_data_pattern = r"---([A-z, ]+)---\n(.+?)\n"

    for i in range(len(inital_mutant_list)):
        with open(f'{data_dir}Mutation_{i}.dat') as f:
            mutant_data_list = re.findall(mutant_data_pattern, f.read())
            for mutant_data in mutant_data_list:
                # prepare muta flag
                muta_flag_list = eval(mutant_data[0])
                muta_flag_list = list(map(lambda x: "".join(x[:1]+x[2:]), muta_flag_list))
                muta_flag = get_unique_mutant(" ".join([initial_mutant_list[i]]+muta_flag_list))            
                # prepare metric data
                metric_data_list = re.findall(metric_data_pattern, mutant_data[1])
                metric_data_mapper = {}
                for metric_data in metric_data_list:
                    metric_data_mapper[metric_data[0]] = eval(metric_data[1])
                data.append((muta_flag, metric_data_mapper))
        break
    return data

data = data_extraction_multiple_initial(initial_mutant_list, proj_path)

#process
for m_flag, m_data in data:
    m_data['BD_norm'] = [i[0] for i in m_data['Bond Dipole']]
    #Unit transfer (e*A0 -> C*cm) & (kcal/(mol*e*A) -> MV/cm) -> * A0/A
    A0 = 5.291*10**-9 #cm
    A = 10**-8 #cm
    #BD
    #BD_norm = BD_norm (* A0 * e)
    BD_norm = np.array(m_data['BD_norm']).astype(float)
    #E
    #E = E * kcal/(A*e*Na*10**6)
    E = np.array(m_data['E']).astype(float)
    m_data['E_mean'] = E.mean()
    #G
    G = -BD_norm * E * A0/A
    m_data['G_mean'] = G.mean()
    del m_data['BD_norm']
    del m_data['Bond Dipole']
    del m_data['E']
    del m_data['traj']
    
# Output
#=========
#---csv---
with open(csv_file,'w') as of:
    # title
    m_name_0, m_data_0 = data[0]
    of.write(f"Mutant,{','.join(m_data_0.keys())}{os.linesep}")
    for m_name, m_data in data:
        of.write(f"{m_name},{','.join(map(lambda x: str(x), m_data.values()))}{os.linesep}")

#---plt---E
    # TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in data]
    # dist = [i['Distance'] for i in data]
    # E_mean = [i['E_mean'] for i in data]
    # BD_mean = [i['BD_norm_mean'] for i in data]
    # p_data = E_mean

    # # E settings
    # plt.xlabel('Distance ()', fontsize=15)
    # plt.ylabel('Field Strength (MV/cm)', fontsize=15)
    # plt.xlim(10,30)

    # plt.scatter(dist, p_data, s=5, c='b')
    # for i, text in enumerate(TAG):
    #     # if text == 'F9A':
    #     #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]-0.3, p_data[i]+0.1))
    #         # continue

    #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]+0.1, p_data[i]+0.1))
    # plt.show()

#---plt---BD
    # TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in data]
    # dist = [i['Distance'] for i in data]
    # E_mean = [i['E_mean'] for i in data]
    # BD_mean = [i['BD_norm_mean'] for i in data]
    # p_data = BD_mean

    # # E settings
    # plt.xlabel('Distance ()', fontsize=15)
    # plt.ylabel('Bond Dipole (a.u.)', fontsize=15)
    # plt.xlim(10,30)

    # plt.scatter(dist, p_data, s=5, c='b')
    # for i, text in enumerate(TAG):
    #     # if text == 'F9A':
    #     #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]-0.3, p_data[i]+0.1))
    #         # continue

    #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]+0.1, p_data[i]))
    # plt.show()

#---plt---G
    # for i in range(len(data)-1,-1,-1):
    #     if data[i]['TAG'] == ('M','A','44','M'):
    #         wt_G_mean = data[i]['G_mean']
    #         del data[i]

    # TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in data]
    # dist = [i['Distance'] for i in data]
    # G_mean = [i['G_mean'] for i in data]
    # p_data = G_mean

    # # E settings
    # plt.xlabel('Distance ()', fontsize=15)
    # plt.ylabel(r'$\bar G$ (kcal/mol)', fontsize=15)
    # plt.xlim(10,30)

    # plt.scatter(dist, p_data, s=5, c='b')
    # for i, text in enumerate(TAG):
    #     # if text == 'F9A':
    #     #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]-0.3, p_data[i]+0.1))
    #         # continue

    #     plt.annotate(text, xy = (dist[i], p_data[i]), xytext = (dist[i]+0.1, p_data[i]+0.1))
    #     plt.plot([10,30], [wt_G_mean, wt_G_mean], c='r')
    # plt.savefig(fname='./Figure_G.svg', format='svg')

#---plt---G-100-hist
    # # TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in data]
    # dist = [i['Distance'] for i in data]
    # G_mean = [i['G_mean'] for i in data]
    # p_data = G_mean

    # # settings
    # plt.xlabel(r'$\bar G$ (kcal/mol)', fontsize=15)
    # plt.ylabel('Frequency', fontsize=15)
    # plt.xlim(-1.5,9)
    # plt.hist(p_data, bins=21, density=0, facecolor='blue', edgecolor='black', alpha=0.7)
    # plt.savefig(fname='./Figure_G_hist.pdf', format='pdf')

    # # probplot(p_data, dist='norm', plot=plt)
    # # plt.savefig(fname='./Figure_G_QQ.png', format='png')

#--xmg-plot---G-100-hist
    # dist = [i['Distance'] for i in data]
    # G_mean = [i['G_mean'] for i in data]
    # p_data = list(G_mean)

    # with open('./G_mean1.dat', 'w') as of:
    #     of.write('#n g_mean\n')
    #     for i, G in enumerate(G_mean):
    #         of.write(str(G)+'\n')

#--pymol-plot--G-100-cartoon
# Index = [i['TAG'][2] for i in data]
# TAG = [i['TAG'] for i in data]

# TAG_map={}

# for i in TAG:
#     TAG_map[i[2]] = None
# for i in TAG:
#     if TAG_map[i[2]] == None:
#         TAG_map[i[2]] = [i[3]]
#     else:
#         TAG_map[i[2]].append(i[3])
#         if TAG_map[i[2]] == i[3]:
#             print('found repeat in: '+i[1]+i[2]+i[3])
# for i in sorted(TAG_map, key=lambda j: int(j)):
#     print(i, sep='', end='+')#,':', TAG_map[i])

#--G metric csv--
# with open('./Mutation_G_100.csv','w') as of:
#     print('TAG', 'G_mean', 'G_SD', 'G_pQ', 'G_mQ', 'G_max', 'G_min', 'G_med', sep=',', file=of)
#     for m_data in data:
#         print(''.join(m_data['TAG']), m_data['G_mean'], m_data['G_SD'], m_data['G_pQ'], m_data['G_mQ'], m_data['G_max'], m_data['G_min'], m_data['G_med'], sep=',', file=of)

#--group by mutation type--
# from AmberMaps import *

# with open('Mutant_group.csv','w') as of:
#     type_list = ['netural', 'charged']
#     # init
#     table9={}
#     for i in type_list:
#         for j in type_list:
#             table9[i+'-'+j] = []
#     # classify
#     for m_data in data:
#         TAG = m_data['TAG']
#         for i in type_list:
#             for j in type_list:
#                 if TAG[0] in resi_subgrp[i] and TAG[3] in resi_subgrp[j]:
#                     table9[i+'-'+j].append(''.join((TAG[0],TAG[2],TAG[3])))
#     # write
#     for i in table9:
#         print(i,','.join(table9[i]),sep=',',end='\n',file=of)
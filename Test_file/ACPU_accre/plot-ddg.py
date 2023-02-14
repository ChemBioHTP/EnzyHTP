
#Set project path here
proj_path='./'
csv_file = f'{proj_path}PuO_1st_ddg.csv'
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

def extract_data_ddg():
    """
    extract ddg data from ddg_result.dat
    """
    result = []
    with open("ddg_results.dat") as f:
        ddg_mapper = eval(f.read().strip())
    for k, v in ddg_mapper.items():
        new_k = []
        existing_r_idx = []
        for k_p in k:
            orig, r_idx, target = k_p.strip().split()
            if r_idx not in existing_r_idx and str(int(r_idx)+resi_idx_mapper["PuO"]["A"]-resi_idx_mapper["PuO"]["B"]) not in existing_r_idx:
                existing_r_idx.append(r_idx)
                new_k.append("".join((orig, r_idx, target)))
        new_k = " ".join(new_k)
        result.append((new_k, v))

    return result

data = extract_data_ddg()

# Output
#=========
#---csv---
with open(csv_file,'w') as of:
    # title
    of.write("Mutant,cart_ddg\n")
    for m_name, m_data in data:
        of.write(f"{m_name},{m_data}\n")

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
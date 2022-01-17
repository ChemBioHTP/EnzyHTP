import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import repeat
from scipy.stats import probplot

# SEt project path here
ProjPath = "./"
Data_file = ProjPath + "Mutation-E-BD-100.dat"
csv_file = ProjPath + "Mutation-E-BD-100.csv"
analysis_list = ["mean", "SD", "max", "min", "med"]
terms = ["BD_norm", "E"]
Data = []

# extract data with eval()
with open(Data_file) as f:
    lines = f.readlines()
    d_flag = 0
    m_flag = 0
    for i, line in enumerate(lines):
        if i == 0:
            m_flag = 1
            m_data = {}
            continue
        if i + 1 == len(lines):
            if d_flag:
                line_data = eval(line.strip())
                m_data[Term] = line_data
                d_flag = 0
            Data.append(m_data)
            break
        if "TAG" in line:
            Data.append(m_data)
            m_flag = 1
            m_data = {}
            continue
        if m_flag:
            MutaFlag = eval(line.strip())[0]
            m_data["TAG"] = MutaFlag
            m_flag = 0
            continue

        if line.strip()[:3] + line.strip()[-3:] == "------":
            Term = line.strip().strip("-")
            d_flag = 1
            continue
        if d_flag:
            line_data = eval(line.strip())
            m_data[Term] = line_data
            d_flag = 0
            continue

# process
for m_data in Data:
    m_data["BD_norm"] = [i[0] for i in m_data["Bond Dipole"]]
    # BD
    BD_norm = np.array(m_data["BD_norm"]).astype(float)
    # Unit transfer (e*A0 -> C*cm)
    A0 = 5.291 * 10 ** -9  # cm
    e = 1.602 * 10 ** -19  # C
    # BD_norm = BD_norm * A0 * e

    m_data["BD_norm_mean"] = BD_norm.mean()
    m_data["BD_norm_SD"] = BD_norm.std()
    m_data["BD_norm_max"] = BD_norm.max()
    m_data["BD_norm_min"] = BD_norm.min()
    m_data["BD_norm_med"] = np.median(BD_norm)

    # E
    E = np.array(m_data["E"]).astype(float)
    # Unit transfer (kcal/(mol*e*A) -> MV/cm)
    A = 10 ** -8  # cm
    Na = 6.02 * 10 ** 23
    kcal = 4184  # J
    # E = E * kcal/(A*e*Na*10**6)

    m_data["E_mean"] = E.mean()
    m_data["E_SD"] = E.std()
    m_data["E_max"] = E.max()
    m_data["E_min"] = E.min()
    m_data["E_med"] = np.median(E)

    # G
    G = (
        -BD_norm * E * A0 / A
    )  # G = -BD_norm * E * A0/A # add minus except for the tests before 2021.7.30.

    m_data["G_mean"] = E.mean()
    m_data["G_SD"] = E.std()
    m_data["G_pQ"] = np.percentile(E, 75)
    m_data["G_mQ"] = np.percentile(E, 25)
    m_data["G_max"] = E.max()
    m_data["G_min"] = E.min()
    m_data["G_med"] = np.median(E)


# Output
# =========
# ---csv---
# with open(csv_file,'w') as of:
#     if_first = 1
#     for m_data in Data:
#         if if_first:
#             line = ''
#             for i in analysis_list:
#                 for j in terms:
#                     line += j+'_'+i+','
#             line+='Distance, TAG\n'
#             of.write(line)
#             if_first = 0
#         line = ''
#         for i in analysis_list:
#             for j in terms:
#                 line += str(m_data[j+'_'+i])+','
#         TAG = m_data['TAG'][0] + m_data['TAG'][2] + m_data['TAG'][3]
#         line+=str(m_data['Distance'])+','+TAG+'\n'
#         of.write(line)

# ---plt---E
# TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in Data]
# dist = [i['Distance'] for i in Data]
# E_mean = [i['E_mean'] for i in Data]
# BD_mean = [i['BD_norm_mean'] for i in Data]
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

# ---plt---BD
# TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in Data]
# dist = [i['Distance'] for i in Data]
# E_mean = [i['E_mean'] for i in Data]
# BD_mean = [i['BD_norm_mean'] for i in Data]
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

# ---plt---G
# for i in range(len(Data)-1,-1,-1):
#     if Data[i]['TAG'] == ('M','A','44','M'):
#         wt_G_mean = Data[i]['G_mean']
#         del Data[i]

# TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in Data]
# dist = [i['Distance'] for i in Data]
# G_mean = [i['G_mean'] for i in Data]
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

# ---plt---G-100-hist
# # TAG = [i['TAG'][0]+i['TAG'][2]+i['TAG'][3] for i in Data]
# dist = [i['Distance'] for i in Data]
# G_mean = [i['G_mean'] for i in Data]
# p_data = G_mean

# # settings
# plt.xlabel(r'$\bar G$ (kcal/mol)', fontsize=15)
# plt.ylabel('Frequency', fontsize=15)
# plt.xlim(-1.5,9)
# plt.hist(p_data, bins=21, density=0, facecolor='blue', edgecolor='black', alpha=0.7)
# plt.savefig(fname='./Figure_G_hist.pdf', format='pdf')

# # probplot(p_data, dist='norm', plot=plt)
# # plt.savefig(fname='./Figure_G_QQ.png', format='png')

# --xmg-plot---G-100-hist
# dist = [i['Distance'] for i in Data]
# G_mean = [i['G_mean'] for i in Data]
# p_data = list(G_mean)

# with open('./G_mean1.dat', 'w') as of:
#     of.write('#n g_mean\n')
#     for i, G in enumerate(G_mean):
#         of.write(str(G)+'\n')

Index = [i["TAG"][2] for i in Data]
TAG = [i["TAG"] for i in Data]

TAG_map = {}

for i in TAG:
    TAG_map[i[2]] = None
for i in TAG:
    if TAG_map[i[2]] == None:
        TAG_map[i[2]] = [i[3]]
    else:
        TAG_map[i[2]].append(i[3])
        if TAG_map[i[2]] == i[3]:
            print("found repeat in: " + i[1] + i[2] + i[3])
for i in sorted(TAG_map, key=lambda j: int(j)):
    print(i, sep="", end="+")  # ,':', TAG_map[i])

# Index.sort(key=lambda i: int(i))
# for i in Index:
#     print(i,sep='', end=',')

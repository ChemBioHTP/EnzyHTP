import collections
from typing import Callable, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import matplotlib
import seaborn as sns
import pickle
import numpy as np

from enzy_htp.mutation_class import Mutation
from enzy_htp import PDBParser

A = 10**-8 #cm
Na = 6.02*10**23
kcal = 4184 #J 
e = 1.602*10**-19 #C     
scale_f = kcal/(A*e*Na*10**6)
font_size_1 = 8
font_size_2 = 10
font_size_3 = 12

def plot_single_stability_violin(
        data_mapper: Dict[Mutation, float], 
        out_path: str, 
        pdb_path: str = None, 
        chain_sync_mapper: Callable=None, 
        merge_chain_sync_mutation: bool=True
    ):
    """make a violin plot for the distribution of single mutations on each position"""
    # prepare data for plot
    site_mapper = collections.defaultdict(list)
    for mut, ddg in data_mapper.items():
        key = mut.get_position_key()
        if chain_sync_mapper and merge_chain_sync_mutation:
            key = chain_sync_mapper(key)
        site_mapper[key].append(ddg)
    if pdb_path:
        stru = PDBParser().get_structure(pdb_path)
        site_mapper = dict(sorted(
            site_mapper.items(), key=lambda x: _get_distance(stru, x[0], ("D", 902))))
        label = []
        for key in site_mapper.keys():
            dist = _get_distance(stru, key, ('D', 902))
            if chain_sync_mapper:
                key = chain_sync_mapper(key, change_chain=False)
            label.append(f"{key[0]}{key[1]}\n({dist:.1f})")
    else:
        site_mapper = dict(sorted(site_mapper.items(), key=lambda x: x[0]))
        label = list(site_mapper.keys())
    x = list(site_mapper.values())

    # make plot
    fig, ax = plt.subplots(figsize=(6.6, 3), dpi=300)
    ax.set_ylim((-15, 30))
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))
    sns.violinplot(data=x, scale="width", width=0.5, linewidth=1, color="orange")
    ax.set_xticks(range(len(label)))
    ax.set_xticklabels(label)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(9)
    ax.set_xlabel('Mutation Site (by distance to substrate)', fontsize=12)
    ax.set_ylabel('Folding Stability (REU)', fontsize=12)
    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def _get_distance(stru, x: Tuple, target: Tuple):
    """get the distance between 2 residues.ca based on the key"""
    p1 = np.array(stru.find_residue_with_key(x).ca_coord)
    p0 = np.array(stru.find_residue_with_key(target).geom_center)
    return np.linalg.norm(p1 - p0)

def plot_stability_vs_ef(mutant_space_data: Dict[Mutation, Tuple[float, float]],
                         out_path: str,
                         s_value: float=0.1, if_add_sele_window: bool=True):
    """takes around 234s"""
    # prepare data
    x = []
    y = []
    for mut, (ef, ddg) in mutant_space_data.items():
        x.append(ef)
        y.append(ddg)

    # make plot
    fig, ax = plt.subplots(figsize=(6.6, 4), dpi=300)
    ax.set_xlim((-10, 140))
    ax.set_ylim((-15, 140))
    ax.scatter(x, y, s=s_value, edgecolor="none")
    if if_add_sele_window:
        add_sele_window(ax, [
            ((min(y), 0), "300"),
            ((0,10), "300"),
            ((10,20), "300"),
            ((20,30), "300"),
            ((30,40), "300"),
            ], x, y)

    # format
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    ax.set_xlabel('Electric Field Strength (MV/cm)', fontsize=font_size_3)
    ax.set_ylabel('Folding Stability (REU)', fontsize=font_size_3)

    fig.tight_layout()
    plt.savefig(out_path, dpi=300)

def add_sele_window(ax: matplotlib.axes.Axes, sele_def: list, data_x, data_y):
    """add selection window to a plot"""
    for sele in sele_def:
        y_low = sele[0][0]
        y_high = sele[0][1]
        if sele[1] == "300":
            xy = zip(data_x, data_y)
            x_sele = map(lambda a: a[0], filter(lambda x: y_low < x[1] <= y_high, xy))
            x_rank = sorted(x_sele, reverse=True)[:300]
            x_low = x_rank[-1]
            x_high = x_rank[0]
        else:
            raise Exception

        width = x_high - x_low
        height = y_high - y_low

        rect = Rectangle(
            (x_low, y_low), width, height,
            linewidth=1, edgecolor='red', facecolor='none')
        ax.add_patch(rect)

# region TODO recycle
def plot_ele_rank(data_path, out_path):
    with open(data_path, "rb") as f:
        ranked_list = pickle.load(f)
    x = [i[1] * scale_f for i in ranked_list]

    fig, ax = plt.subplots(figsize=(3.3, 2), dpi=300)
    ax.set_xlim((0, 150))
    ax.set_ylim((0, 0.02))

    # ax.hist(x, bins=30, alpha=0.5, color='blue')
    sns.kdeplot(x=x, ax=ax,
                color='orange', fill=True, alpha=0.5)

    # 设置坐标轴上数字的小数点位数
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('Electric Field Strength (MV/cm)', fontsize=font_size_1)
    ax.set_ylabel('Probability Density', fontsize=font_size_1)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_stability_rank(data_path, out_path):
    with open(data_path, "rb") as f:
        ranked_list = pickle.load(f)
    x = [i[1] for i in ranked_list]

    fig, ax = plt.subplots(figsize=(3.3, 2), dpi=300)
    # ax.set_xlim((0, 150))
    # ax.set_ylim((0, 0.02))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))

    # ax.hist(x, bins=30, alpha=0.5, color='blue')
    sns.kdeplot(x=x, ax=ax,
                color='orange', fill=True, alpha=0.5)

    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('Folding Stability (REU)', fontsize=font_size_1)
    ax.set_ylabel('Probability Density', fontsize=font_size_1)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_stability_vs_ef_good_ddg(data_path, out_path, s_value=0.1):
    with open(data_path, "rb") as f:
        ef_stability_rank_list = pickle.load(f)
    x = []
    y = []
    for mut, ef, ddg in ef_stability_rank_list:
        x.append(ef*scale_f)
        y.append(ddg)

    fig, ax = plt.subplots(figsize=(3.3, 2), dpi=300)
    # ax.set_xlim((0, 150))
    # ax.set_ylim((0, 0.02))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))

    ax.scatter(x, y, s=s_value, edgecolor="none")
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('Electric Field Strength (MV/cm)', fontsize=font_size_1)
    ax.set_ylabel('Folding Stability (REU)', fontsize=font_size_1)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_ef_vs_num_mutation(ef_data, out_path, s_value=0.1):
    """scatter plot of ef change vs number of mutations
    took 135.32891771s"""
    with open(ef_data, "rb") as f:
        ef_rank_list = pickle.load(f)
    x = []
    y = []
    for mut, ef in ef_rank_list:
        x.append(ef*scale_f)
        y.append(len(mut))

    fig, ax = plt.subplots(figsize=(6.6, 4), dpi=300)
    ax.set_xlim((0, 140))
    # ax.set_ylim((0, 0.02))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

    ax.scatter(x, y, s=s_value, edgecolor="none")
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('Electric Field Strength (MV/cm)', fontsize=font_size_2)
    ax.set_ylabel('Number of mutations', fontsize=font_size_2)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_ef_vs_num_mutation_violin(ef_data, out_path):
    with open(ef_data, "rb") as f:
        ef_rank_list = pickle.load(f)
    x = range(13)
    y = [[] for _ in range(13)]
    for mut, ef in ef_rank_list:
        y[len(mut)].append(ef*scale_f)
    fig, ax = plt.subplots(figsize=(6.6, 4), dpi=300)
    # ax.set_xlim((0, 140))
    # ax.set_ylim((0, 0.02))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(10))

    sns.violinplot(data=y, scale="width", width=0.5, linewidth=1, color="orange")
    # ax.set_xticks(range(len(x)))
    # ax.set_xticklabels(x)
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('Number of Mutations', fontsize=font_size_2)
    ax.set_ylabel('Electric Field Strength (MV/cm)', fontsize=font_size_2)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_ef_vs_num_mutation_violin_good_ddg(ef_data, out_path):
    with open(ef_data, "rb") as f:
        ef_ddg_rank_list = pickle.load(f)
    x = range(13)
    y = [[] for _ in range(13)]
    for mut, ef, ddg in ef_ddg_rank_list:
        y[len(mut)].append(ef*scale_f)
    print(y[2])
    fig, ax = plt.subplots(figsize=(6.6, 4), dpi=300)
    # ax.set_xlim((0, 140))
    ax.set_ylim((0, 110))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))

    sns.violinplot(data=y, scale="width", width=0.5, linewidth=1, color="orange")
    # ax.set_xticks(range(len(x)))
    # ax.set_xticklabels(x)
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    ax.set_xlabel('Number of Mutations', fontsize=12)
    ax.set_ylabel('Electric Field Strength (MV/cm)', fontsize=12)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

def plot_verify_sample_points(sample_data_list, stablilty_data, ef_data, out_path, s_value=10, color="red"):
    """the sample points plotted as aligned with the stability_vs_ef plot"""
    with open(ef_data, "rb") as f:
        ef_rank_list = pickle.load(f)
    with open(stablilty_data, "rb") as f:
        stablilty_rank_list = pickle.load(f)
    stablilty_rank_mapper = dict(stablilty_rank_list)
    ef_rank_mapper = dict(ef_rank_list)
    x = []
    y = []
    for mut, value in sample_data_list:
        x.append(ef_rank_mapper[mut]*scale_f)
        y.append(stablilty_rank_mapper[mut])

    fig, ax = plt.subplots(figsize=(6.6, 4), dpi=300)
    ax.set_xlim((-10, 140))
    ax.set_ylim((-15, 140))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.005))

    ax.scatter(x, y, s=s_value, edgecolor="none", c=color)
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_2)
    ax.set_xlabel('Electric Field Strength (MV/cm)', fontsize=font_size_3)
    ax.set_ylabel('Folding Stability (REU)', fontsize=font_size_3)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300, transparent=True)


def _make_good_stability_ef_csv(data_path, out_path):
    import csv
    from enzy_htp.mutation.mutation import get_mutant_name_tag
    with open(data_path, "rb") as f:
        result = pickle.load(f)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mutant", "EF", "ddG"])
        writer.writerows([[get_mutant_name_tag(i)[1:].replace("_"," "),j*scale_f,k] for i,j,k in result])

def plot_2nd_metric_dist(data_1: list, data_2: list, out_path, s_value=0.1, if_add_sele_window=0):
    """plot the 2d metric distribution"""
    x = data_1
    y = data_2

    fig, ax = plt.subplots(figsize=(3.3, 2), dpi=300)
    ax.set_xlim((-2.3, 1.2))
    ax.set_ylim((0.7,2.2))
    # 设置坐标轴上数字的小数点位数
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    # 设置 X 和 Y 轴上的刻度数量
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))

    ax.scatter(x, y, s=s_value, edgecolor="none")
    if if_add_sele_window:
        add_sele_window(ax, [
            ((min(y), 0), "300"),
            ((0,10), "300"),
            ((10,20), "300"),
            ((20,30), "300"),
            ((30,40), "300"),
            ], x, y)
    ax.axhline(WT_METRIC["sasa_avg"]-0.2, color='grey', linestyle='dashed', linewidth=1)
    ax.axhline(WT_METRIC["sasa_avg"]+0.3, color='grey', linestyle='dashed', linewidth=1)
    # # 获取 X 和 Y 轴上的 Tick 对象，并设置数字的字体大小
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(font_size_1)
    ax.set_xlabel('$G_{elec,C-H}$ (kcal/mol)', fontsize=font_size_2)
    # ax.set_ylabel('Folding Stability (REU)', fontsize=font_size_2)
    # ax.set_ylabel('pKa', fontsize=font_size_2)
    ax.set_ylabel('SPI', fontsize=font_size_2)

    fig.tight_layout() 
    plt.savefig(out_path, dpi=300)

# endregion

def puo_chain_sync(key: Tuple[str, int], change_chain: bool=True) -> Tuple[str, int]:
    """desc dimer in PuO"""
    if key[0] == "B":
        new_idx = key[1] - 450
        if change_chain:
            new_chain = "A"
        else:
            new_chain = "B"
        return (new_chain, new_idx)
    else:
        return key

def main():
    with open("mutant_space_1_coarse_def_ddg_fold.pickle", "rb") as f:
        mutant_space_data = pickle.load(f)
    with open("single_mutation_ddg.pickle", "rb") as f:
        single_mutation_data = pickle.load(f)
    plot_single_stability_violin(single_mutation_data, "figures/single_stability_violin.svg", pdb_path="puo_acp.pdb", 
                                 chain_sync_mapper=puo_chain_sync, merge_chain_sync_mutation=False)
    plot_single_stability_violin(single_mutation_data, "figures/single_stability_violin_sync.svg", pdb_path="puo_acp.pdb",
                                 chain_sync_mapper=puo_chain_sync)
    plot_stability_vs_ef(mutant_space_data, "figures/mut_space_1_ddgfold_vs_ef.svg", s_value=0.1, if_add_sele_window=False)

    pass

if __name__ == "__main__":
    main()
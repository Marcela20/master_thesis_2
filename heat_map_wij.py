"""
Script to generate heatmap with RSCU values.

"""

import seaborn as sns
from CAI import RSCU
import pandas as pd
import matplotlib.pyplot as plt
from utils import create_list_of_records
import os

codons = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0, "CTT": 0,
    "CTC": 0, "CTA": 0, "CTG": 0, "ATT": 0, "ATC": 0,
    "ATA": 0, "ATG": 0, "GTT": 0, "GTC": 0, "GTA": 0,
    "GTG": 0, "TAT": 0, "TAC": 0, "CAT": 0, "CAC": 0,
    "CAA": 0, "CAG": 0, "AAT": 0, "AAC": 0, "AAA": 0,
    "AAG": 0, "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0, "CCT": 0,
    "CCC": 0, "CCA": 0, "CCG": 0, "ACT": 0, "ACC": 0,
    "ACA": 0, "ACG": 0, "GCT": 0, "GCC": 0, "GCA": 0,
    "GCG": 0, "TGT": 0, "TGC": 0, "TGG": 0, "CGT": 0,
    "CGC": 0, "CGA": 0, "CGG": 0, "AGT": 0, "AGC": 0,
    "AGA": 0, "AGG": 0, "GGT": 0, "GGC": 0, "GGA": 0,
    "GGG": 0}
synonymousCodons = {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"]
}


def create_list_of_codons(dict_of_codons):
    list_codons = []
    [list_codons.append(k) for k in dict_of_codons]
    list_codons = [codon for codon in dict_of_codons]
    return list_codons


def create_dict_of_dicts():
    """
    returns nested dict like {"TRP":["TGG":0...]...} where codons have 0 value
    """
    for key, value in synonymousCodons.items():
        value2 = dict.fromkeys(value, 0)
        synonymousCodons[key] = value2
    return synonymousCodons


def order_codons(codons_list):
    c_ending = []
    g_ending = []
    a_ending = []
    t_ending = []

    for i in codons_list:

        if i[-1] == "C":
            c_ending.append(i)
        if i[-1] == "G":
            g_ending.append(i)
        if i[-1] == "A":
            a_ending.append(i)
        if i[-1] == "T":
            t_ending.append(i)
    list_of_all = c_ending[::-1] + g_ending[::-1] + a_ending[::-1] + t_ending[::-1]
    return list_of_all


def get_RSCU(list_of_dicts):
    list_of_all = order_codons(create_list_of_codons(codons))
    seqs = []

    for gene in list_of_dicts:
        seqs.append(gene['seq'])
    RSCU_ref_data_not_ordered = RSCU(seqs)  # dictionary
    RSCU_ref_data = {k: RSCU_ref_data_not_ordered[k] for k in list_of_all if k not in ["ATG", "TGG"]}

    return [RSCU_ref_data]


def get_wij(dict_of_rscu):
    list_of_all = order_codons(create_list_of_codons(codons))
    dict_of_wij_vals = {}
    sorted_dict_of_wij_vals = {}
    codons_for_genom = create_dict_of_dicts()  # {"TRP":["TGG":0...]...}
    for aminoacid, dict_codons in codons_for_genom.items():
        for codon, value in dict_of_rscu.items():
            if codon in dict_codons:
                dict_codons[codon] = dict_of_rscu[codon]
        RSCUimax = max(dict_codons.values())
        if aminoacid in ['MET', 'TRP']:
            continue
        for codon, rscu_value in dict_codons.items():
            dict_of_wij_vals[codon] = rscu_value / RSCUimax

    RSCU_ref_data = {k: dict_of_wij_vals[k] for k in list_of_all if k not in ["ATG", "TGG"]}

    return RSCU_ref_data


def generate_tuple_for_height_ratios(num_of_groups):
    start_list = []
    for i in range(len(num_of_groups)):
        start_list.append(.9)
    start_list.append(.3)
    start_list = tuple(start_list)
    return start_list


def get_mean_RSCU_of_many_files(directory):
    counter = 0
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if counter > 0:
            df_single = pd.DataFrame(get_RSCU(create_list_of_records(f'{directory}\{filename}', 200)))
            # df_single = pd.DataFrame(dict(zip(df_single.columns, list(df_single.mean()))), index=[0])
            df = df.append(df_single, ignore_index=True)


        else:
            # print(get_RSCU(create_list_of_records(f'{directory}\{filename}', 200)))
            df = pd.DataFrame(get_RSCU(create_list_of_records(f'{directory}\{filename}', 200)))
            df = pd.DataFrame(dict(zip(df.columns, list(df.mean()))), index=[0])

        counter += 1
    return df


if __name__ == "__main__":
    directiories = ['GC_high', 'GC_medium', 'GC_low']
    data_frames_to_plot = []
    tup_of_height_ratios = generate_tuple_for_height_ratios(directiories)
    num_of_subplots = len(directiories) + 1
    for directory in directiories:
        data_frames_to_plot.append(get_mean_RSCU_of_many_files(f'K:\{directory}'))

    grid_kws = {"height_ratios": tup_of_height_ratios, "hspace": .3}

    fig, (ax1, ax2, ax7, ccbar) = plt.subplots(num_of_subplots, figsize=(13, 8), gridspec_kw=grid_kws,
                                               constrained_layout=True)  # add here before ax7, ax(1, 2, 3..) dependinding on number of organisms, if You have one, ax7 is sufficient
    l_of_axes = [ax1,
                 ax2,
                 ax7]  # add here before ax7, ax(1, 2, 3..) dependinding on number of organisms, if You have one, ax7 is sufficient
    counter = 0
    for df, axe in zip(data_frames_to_plot, l_of_axes):
        counter += 1
        if axe == ax7:
            heat = sns.heatmap(df, yticklabels=False, xticklabels=True, vmin=0, vmax=6, cmap="gist_heat", ax=axe,
                               cbar_ax=ccbar,
                               cbar_kws={"orientation": "horizontal"}, linewidths=0)
        else:
            heat = sns.heatmap(df, yticklabels=False, xticklabels=False, cmap="gist_heat", ax=axe,
                               cbar_ax=ccbar, cbar_kws={"orientation": "horizontal"}, linewidths=0)
        file_name = 'GC content'
        if counter == 0:
            heat.set_ylabel(f'high {file_name}')
        elif counter == 1:
            heat.set_ylabel(f'medium {file_name}')
        else:
            heat.set_ylabel(f'low {file_name}')

    plt.show()

    # save_heatmaps_to_folder(list_of_files, f' ')

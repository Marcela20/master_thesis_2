import pandas as pd
import matplotlib.pyplot as plt
from ramps_cai import get_fixed_size_windows_from_seq
from utils import create_list_of_records


def get_absolute_adaptiveness_values(file):
    dict_1 = {}
    with open(file, 'r') as wi_values:
        a = wi_values.readlines()
    res = map(lambda x: x.split(), a)  # nested list

    def dict_of_codons(x, y):
        return dict_1.setdefault(x, y)

    for inner_list in res:
        if len(inner_list) == 2:
            abc = inner_list[1].replace(",", ".")
            dict_of_codons(inner_list[0], abc)
        else:
            dict_of_codons(inner_list[1], inner_list[2])

    return dict_1



def check_wij_for_ramps(list_of_ramps, wij_of_codons):
    """returns sum of all wij values for each ramp for one gene"""

    ramps = []  # lista sumy wij dla kodonów we wszystkich rampach
    if list_of_ramps is not None:
        for ramp in list_of_ramps:

            val_ramp = 0
            for i in range(0, len(ramp), 3):
                val_ramp += float(wij_of_codons[ramp[i:i + 3]])
            ramps.append(val_ramp / (len(ramp) / 3))

        return ramps
    else:
        print("this was None")


def check_wij_for_ramps_for_separate_plots(list_of_ramps, wij_of_codons):
    """returns ramps of only these genes where there is a ramp of poorly adapted codons at the beginning of a gene"""
    # dict like {codon: wij value}
    ramps = []  # lista sumy wij dla kodonów we wszystkich rampach

    for ramp in list_of_ramps:

        val_ramp = 0
        for i in range(0, len(ramp), 3):
            val_ramp += float(wij_of_codons[ramp[i:i + 3]])
        ramps.append(val_ramp / (len(ramp) / 3))
    if sum(ramps[0:len(ramps) // 2]) < sum(ramps[len(ramps) // 2:-1]):
        return ramps


def create_plots_of_all_genes(sequences, path, codons):
    """creates plots of only these genes that have a ramps of poorly adapted codons"""
    for seq in sequences:
        ramps = get_fixed_size_windows_from_seq(seq, len_of_ramps, step, num_of_ramps)
        wij_of_ramps = check_wij_for_ramps_for_separate_plots(ramps, codons)

        if wij_of_ramps != None:
            df = pd.DataFrame(wij_of_ramps)
            file = str(seq["name"]).replace("|", "_")

            plt.plot(df)
            plt.title(seq["protein"])
            plt.savefig(path + f"{file}.png")
            plt.close()


def create_plot_of_sum_of_genes(sequences, path, name, number_of_ramps, codons):
    """change value of number_of_ramps when You change len of wind or step in function
    'get_fixed_size_windows_from_seq' return plot of sums of RSCUij for ramps in whole genom (all genes, with poorly
    adapted codons at the beginning and without) """

    dummy_data = {0: [0] * number_of_ramps}
    ultimate_data_frame = pd.DataFrame(dummy_data)

    for seq in sequences:
        ramps = get_fixed_size_windows_from_seq(seq, len_of_ramps, step, num_of_ramps)
        wij_of_ramps = check_wij_for_ramps(ramps, codons)
        df = pd.DataFrame(wij_of_ramps)
        ultimate_data_frame = ultimate_data_frame.add(df, fill_value=0)

    ultimate_data_frame = ultimate_data_frame.div(len(sequences))
    file = name.replace(".txt", "")
    plt.plot(ultimate_data_frame)
    plt.ylabel("wij value")
    plt.xlabel("ramp")
    plt.savefig(path + f"{file}_{len_of_ramps}_tuller_tAI.png")
    plt.close()


if __name__ == '__main__':
    list_of_files = ["sample_cds_fasta.txt"]
    path_of_folder = "C:\\Users\\marce\\Desktop\\ramps_of_all_sum_2\\"

    len_of_ramps = 12
    step = 3
    num_of_ramps = 200
    len_of_gene = 3 + len_of_ramps + (step * (num_of_ramps - 1)) + 1
    codons_wij_vals = get_absolute_adaptiveness_values(
        'old_tAIs/idk.txt')  # this file must be adjusted to right organism
    print(len(codons_wij_vals))



    for file in list_of_files[0:1]:
        organism = f"sequences_of_model_organisms/{file}"
        sequences_of_an_organism = create_list_of_records(organism, len_of_gene)
        create_plot_of_sum_of_genes(sequences_of_an_organism, path_of_folder, file, num_of_ramps, codons_wij_vals)

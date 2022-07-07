from Bio.SeqUtils import GC123
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from utils import create_list_of_records
from ramps_cai import get_fixed_size_windows_from_seq
import os


def count_GC_in_each_position(list_of_dicts):
    initial_dataframe = pd.DataFrame([(0, 0, 0, 0)], columns=['total', 'first_pos', 'second_pos', 'third_pos'])
    for dict_for_gene in list_of_dicts:
        df = pd.DataFrame([GC123(dict_for_gene["seq"])], columns=['total', 'first_pos', 'second_pos', 'third_pos'])
        initial_dataframe = initial_dataframe.add(df, fill_value=0)
    initial_dataframe = initial_dataframe / len(list_of_dicts)
    return initial_dataframe


def count_GC_in_ramps(list_of_dicts):
    list_of_GC_values = []
    df = None
    for gene in list_of_dicts:
        list_of_GC_values = []
        list_of_ramps = get_fixed_size_windows_from_seq(gene, len_of_ramps, step, num_of_ramps)

        for ramp in list_of_ramps:
            list_of_GC_values.append(list(GC123(ramp)))

        if df is None:
            df = pd.DataFrame(list_of_GC_values, columns=['total', 'first pos', 'second pos', 'third pos'])

        else:
            df_new = pd.DataFrame(list_of_GC_values, columns=['total', 'first pos', 'second pos', 'third pos'])
            df = df.add(df_new, fill_value=0)

    df = df / len(list_of_dicts)
    return df


def count_gc_for_many_files(path_to_folder, path_dest_folder, name_of_file, path_for_overall=None, path_for_ramps=None):
    directory = os.fsencode(path_to_folder)
    dummy_data = {'total': 0, 'first_pos': 0, 'second_pos': 0, 'third_pos': 0}
    df_overall_GC = pd.DataFrame(dummy_data, index=[0])

    counter = 0
    for file in os.listdir(directory):

        filename = os.fsdecode(file)

        organism = f"{directory.decode('ascii')}\{filename}"

        list_of_dicts_of_cds = create_list_of_records(organism, len_of_gene)
        df_overall_GC = df_overall_GC.add(count_GC_in_each_position(list_of_dicts_of_cds), fill_value=0)
        if counter == 0:
            data_to_plot = count_GC_in_ramps(list_of_dicts_of_cds).div(100000000)
        else:
            subdata_to_plot = count_GC_in_ramps(list_of_dicts_of_cds).div(100000000)
            data_to_plot = data_to_plot.add(subdata_to_plot)

        counter += 1
    df_overall_GC = df_overall_GC.div(counter)
    data_to_plot = data_to_plot.div(counter)
    data_to_plot = data_to_plot.multiply(100000000)

    print(df_overall_GC)
    df_overall_GC.to_csv(path_for_overall)
    data_to_plot.to_csv(path_for_ramps)
    plot = data_to_plot.plot()
    plot.set_ylabel("percentage of GC (%)")
    plot.set_xlabel("ramp")
    plt.legend(loc="upper center", bbox_to_anchor=(0.74, 1.15), ncol=2)
    locs, labels = plt.xticks()
    plt.xticks(np.arange(0, num_of_ramps + 1, step=20))
    plt.yticks(np.arange(0, 100, step=10))

    plt.rc('xtick', labelsize=30)
    plt.rc('ytick', labelsize=30)
    plt.savefig(path_dest_folder + f"ultimate_sum_GC_{name_of_file}.png")
    plt.close()


if __name__ == "__main__":
    len_of_ramps = 3
    step = 3
    num_of_ramps = 200
    len_of_gene = 3 + len_of_ramps + (step * (num_of_ramps - 1)) + 1

    count_gc_for_many_files(r'K:\GC_low', r"C:\\Users\\marce\\Desktop\\GC_percent\\", "low",
                            r'to_appendix/GC_low_overall', r'to_appendix/GC_low_ramps')
from CAI import RSCU
import pandas as pd
import matplotlib.pyplot as plt
from utils import performance, create_list_of_records_no_ref_set
import os

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


def create_dict_of_dicts():
    """
    returns nested dict like {"TRP":["TGG":0...]...} where codons have 0 value
    """
    for key, value in synonymousCodons.items():
        value2 = dict.fromkeys(value, 0)
        synonymousCodons[key] = value2
    return synonymousCodons


def get_RSCU(list_of_sequences):
    """
    returns dict like {codon:value...}
    """
    seqs = []
    for dictionary in list_of_sequences:
        seqs.append(dictionary["seq"])

    RSCU_values = RSCU(seqs, genetic_code=11)
    return RSCU_values


def get_RSCUij(dictionary_of_RSCU_vals, dictionary_of_RSCU_vals_reference_set):
    """
    returns nested dictionary {aminoacid:{codon:wij value, codon:wij value}....}
    each codon has value of it's RSCU divided by RSCU of most common one from a reference set
    """

    codons_for_genom = create_dict_of_dicts()

    dictionary_of_RSCU_vals = dict(dictionary_of_RSCU_vals)  # {codon:RSCU_value, ...}
    dictionary_of_RSCU_vals_reference_set = dict(dictionary_of_RSCU_vals_reference_set)  # {"TRP":["TGG":RSCU...]...}

    for key, value in codons_for_genom.items():
        for key_2, val_2 in dictionary_of_RSCU_vals.items():
            if key_2 in value:
                value[key_2] = dictionary_of_RSCU_vals[key_2]

    for aa, v in dictionary_of_RSCU_vals_reference_set.items():  # v is dict of codons
        max_val = max(v.values())
        for codon, RSCU_val in codons_for_genom[aa].items():
            codons_for_genom[aa][codon] = RSCU_val / max_val

    return codons_for_genom


def get_RSCUij_for_ref_set(dictionary_of_RSCU_vals, dict_of_codons):
    """
    returns nested dictionary {aminoacid:{codon:wij value, codon:wij value}....}
    each codon has value of it's RSCu devided by RSCu of most common one
    """
    codons = dict_of_codons.copy()
    dictionary_of_RSCU_values = dict(dictionary_of_RSCU_vals)
    for key, value in codons.items():
        for key_2, val_2 in dictionary_of_RSCU_values.items():
            if key_2 in value:
                value[key_2] = dictionary_of_RSCU_values[key_2]

    return codons


def get_fixed_size_windows_from_seq(sequence, len_of_wind, step_size, num_of_occurrences):
    """
    returns ramps for one gene
    """
    if step_size % 3 != 0:
        raise Exception('Step must be multiply of 3')
    if len_of_wind % 3 != 0:
        raise Exception('len_of_wind must be multiply of 3')
    occurrences = []
    counter = 0
    for num_nt in range(3, len(sequence["seq"]) - len_of_wind, step_size):
        if counter <= num_of_occurrences:
            counter += 1
            occurrence = sequence["seq"][num_nt: num_nt + len_of_wind]
            occurrences.append(occurrence)
    return occurrences



def check_wij_for_ramps(list_of_ramps, wij_of_codons):
    """
    returns sum of all RSCUij values for each ramp for one gene
    """
    codons_2 = {}
    ramps = []  # lista sumy wij dla kodonów we wszystkich rampach
    for key, val in wij_of_codons.items():
        for k, v in val.items():
            codons_2[k] = v
    for ramp in list_of_ramps:
        val_ramp = 0
        for i in range(0, len(ramp), 3):
            val_ramp += codons_2[ramp[i:i + 3]]
        ramps.append(val_ramp / (len(ramp) / 3))
    return ramps


def check_wij_for_ramps_for_separate_plots(list_of_ramps, wij_of_codons):
    """
    returns ramps of only these genes where there is a ramp of poorly adapted codons at the beginning of a gene
    """
    codons_2 = {}  # dict like {codon: wij value}
    ramps = []  # lista sumy wij dla kodonów we wszystkich rampach
    for key, val in wij_of_codons.items():
        for k, v in val.items():
            codons_2[k] = v
    for ramp in list_of_ramps:
        val_ramp = 0
        for i in range(0, len(ramp), 3):
            val_ramp += codons_2[ramp[i:i + 3]]
        ramps.append(val_ramp / (len(ramp) / 3))
    if sum(ramps[0:len(ramps) // 2]) < sum(ramps[len(ramps) // 2:-1]):
        return ramps


def create_plots_of_all_genes(sequences, path, codons, len_of_ramps, step, num_of_ramps):
    """
    creates plots of only these genes that have a ramps of poorly adapted codons
    """
    for seq in sequences:
        ramps = get_fixed_size_windows_from_seq(seq, len_of_ramps, step, num_of_ramps)
        wij_of_ramps = check_wij_for_ramps_for_separate_plots(ramps, codons)
        if wij_of_ramps is not None:
            df = pd.DataFrame(wij_of_ramps)
            file_name = str(seq["name"]).replace("|", "_")
            plt.plot(df)
            plt.title(seq["protein"])
            plt.savefig(path + f"{file_name}.png")
            plt.close()


def create_plot_of_sum_of_genes(sequences, path, name, number_of_ramps, codons, len_of_ramps, step):
    """
    change value of number_of_ramps when You change len of wind or step in function
    'get_fixed_size_windows_from_seq' return plot of sums of RSCUij for ramps in whole genom (all genes, with poorly
    adapted codons at the beginning and without)
    """
    dummy_data = {0: [0] * number_of_ramps}
    ultimate_data_frame = pd.DataFrame(dummy_data)
    for seq in sequences:
        ramps = get_fixed_size_windows_from_seq(seq, len_of_ramps, step, number_of_ramps)
        wij_of_ramps = check_wij_for_ramps(ramps, codons)
        df = pd.DataFrame(wij_of_ramps)
        ultimate_data_frame = ultimate_data_frame.add(df, fill_value=0)
    ultimate_data_frame = ultimate_data_frame.div(len(sequences))
    file_name = name.replace(".txt", "")
    plt.plot(ultimate_data_frame)
    plt.savefig(path + f"{file_name}_{len_of_ramps}_{number_of_ramps}.png")
    plt.close()

@performance
def create_sum_of_files_cai(num_of_files, directory, dict_of_codons, path_to_export=None):
    path_of_folder = r"C:\\Users\\marce\\Desktop\\ramps_of_all_sum_2\\"
    len_of_ramps = 3
    step = 3
    num_of_ramps = 200
    len_of_gene = 3 + len_of_ramps + (step * (num_of_ramps - 1)) + 1
    directory = os.fsencode(directory)
    dummy_data = {0: [0] * num_of_ramps}
    ultimate_data_frame = pd.DataFrame(dummy_data)
    dummy_data_for_sum = {0: [0] * num_of_files}
    ultimate_data_frame_sum = pd.DataFrame(dummy_data_for_sum)

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        organism = f"{directory.decode('ascii')}\{filename}"
        sequences_of_an_organism_long = []
        sequences_of_an_organism = create_list_of_records_no_ref_set(organism, 3)

        for dict_seq in sequences_of_an_organism[0]:
            if len(dict_seq['seq']) >= len_of_gene:
                sequences_of_an_organism_long.append(dict_seq)

        seqs_ref_set = sequences_of_an_organism[1]
        RSCU_of_ref_set = get_RSCUij_for_ref_set(get_RSCU(seqs_ref_set), dict_of_codons)
        codons_RSCU = get_RSCUij(get_RSCU(sequences_of_an_organism_long), RSCU_of_ref_set)

        for seq in sequences_of_an_organism_long:
            ramps = get_fixed_size_windows_from_seq(seq, len_of_ramps, step, num_of_ramps)
            wij_of_ramps = check_wij_for_ramps(ramps, codons_RSCU)
            df = pd.DataFrame(wij_of_ramps)
            ultimate_data_frame = ultimate_data_frame.add(df, fill_value=0)
        ultimate_data_frame = ultimate_data_frame.div(len(sequences_of_an_organism_long))
        ultimate_data_frame_sum = ultimate_data_frame_sum.add(ultimate_data_frame, fill_value=0)
    ultimate_data_frame_sum = ultimate_data_frame_sum.div(num_of_files)
    ultimate_data_frame.to_csv(path_to_export)
    plt.plot(ultimate_data_frame_sum)
    plt.savefig(path_of_folder + f"ultimate_sum_{len_of_ramps}_{num_of_ramps}low.png")
    plt.close()


if __name__ == "__main__":
    DICT_OF_CODONS = create_dict_of_dicts()
    create_sum_of_files_cai(101, r'K:\GC_low', r'to_appendix/GC_low', DICT_OF_CODONS)

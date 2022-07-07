from Bio import SeqIO
import re
from utils import create_list_of_records


def get__ribosomal_protein_from_description(ref_list_of_dicts):
    reference_set = []
    for dictionary in ref_list_of_dicts:
        if_ribosomal = re.search('ribosomal', dictionary["protein"])
        if if_ribosomal:
            if_mitochondrial = re.search('mitochondrial', dictionary["protein"])
            if not if_mitochondrial:
                # if_subunit = re.search('subunit', dictionary["protein"])
                # if if_subunit:
                reference_set.append(dictionary)
    print('ribosomal proteins',reference_set)
    return reference_set



if __name__ == "__main__":

    get__ribosomal_protein_from_description(create_list_of_records("sequences_of_model_organisms/S. cerevisiae.txt", 304))

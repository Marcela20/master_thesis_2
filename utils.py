from Bio import SeqIO
import re
from time import time


def performance(fn):
    def wrapper(*args, **kwargs):
        time_0 = time()
        result = fn(*args, **kwargs)
        time_1 = time()
        perf = time_1 - time_0
        print(f'{fn} took {perf} sek')
        return result

    return wrapper


@performance
def create_list_of_records(file, min_len_of_gene):
    ref_list = []
    stop_codons = ["TAG", "TGA", "TAA"]
    record_ref = SeqIO.parse(file, "fasta")

    def is_stop_codon(x):
        return any(x[i:i + 3] in stop_codons for i in range(0, len(x), 3))

    for r in record_ref:
        r_len = len(r.seq)
        if r_len >= min_len_of_gene and r_len % 3 == 0:

            try:
                protein_name = re.search(r'\[protein=(.*?)\]', r.description).group(1)
                if protein_name == "hypothetical protein":
                    continue

            except AttributeError:
                continue
            if r.seq[-3:] in stop_codons:
                if re.match('^[ATCG]+$', str(r.seq)) and is_stop_codon(r.seq[0:-3]) is False:
                    ref_list.append({
                        "seq": str(r.seq),
                        "protein": protein_name,
                        "name": r.name})

    return ref_list

@performance
def create_list_of_records_no_ref_set(file, min_len_of_gene):
    ref_list = []
    reference_genes = []
    stop_codons = ["TAG", "TGA", "TAA"]
    record_ref = SeqIO.parse(file, "fasta")

    def is_stop_codon(x):
        return any(x[i:i + 3] in stop_codons for i in range(0, len(x), 3))

    for r in record_ref:
        r_len = len(r.seq)
        if r_len >= min_len_of_gene and r_len % 3 == 0:
            try:
                protein_name = re.search(r'\[protein=(.*?)\]', r.description).group(1)
                if protein_name == "hypothetical protein":
                    continue
                if_ribosomal = re.search('ribosomal', protein_name)
                if if_ribosomal:
                    if_mitochondrial = re.search('mitochondrial', protein_name)
                    if not if_mitochondrial:
                        if r.seq[-3:] in stop_codons:
                            if re.match('^[ATCG]+$', str(r.seq)) and is_stop_codon(r.seq[0:-3]) is False:
                                reference_genes.append({
                                    "seq": str(r.seq),
                                    "protein": protein_name,
                                    "name": r.name})

            except AttributeError:
                continue
            if r.seq[-3:] in stop_codons:
                if re.match('^[ATCG]+$', str(r.seq)) and is_stop_codon(r.seq[0:-3]) is False:
                    ref_list.append({
                        "seq": str(r.seq),
                        "protein": protein_name,
                        "name": r.name})

    return [ref_list, reference_genes]
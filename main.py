"""generating ramps with GC percent, RSCU heatmap and ramps with CAI values"""

from RSCU_for_many_files import *
from ramps_cai import *
from ramps_tAI import *
from heat_map_wij import *
from read_csv import *
from get_ref_set import *
from percent_GC import *


if __name__ == '__main__':

    len_of_ramps = 3
    step = 3
    num_of_ramps = 200
    len_of_gene = 3 + len_of_ramps + (step * (num_of_ramps - 1)) + 1
    directiories = ['GC_high', 'GC_medium', 'GC_low']
    DICT_OF_CODONS = create_dict_of_dicts()

    count_gc_for_many_files(br"GC_organisms\\GC_low", r"appendix\\GC_percent\\", "low",  len_of_gene, num_of_ramps, len_of_ramps, step,
                            r'to_appendix/GC_low_overall.csv', r'to_appendix/GC_low_ramps.csv')

    create_RSCU_heatmap(directiories, "GC_organisms/")

    create_sum_of_files_cai(101, r'GC_organisms\\GC_high', DICT_OF_CODONS, r'to_appendix/GC_high.csv')

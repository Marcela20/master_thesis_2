from utils import *
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

    count_gc_for_many_files(br"C:\\Users\\marce\Desktop\\master_thesis_project\\master_thesis_2_0\\GC_organisms\\GC_high", r"C:\\Users\\marce\\Desktop\\GC_percent\\", "low",
                            r'to_appendix/GC_low_overall', r'to_appendix/GC_low_ramps')

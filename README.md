### Scripts for my master thesis project regarding Codon Usage Bias (CUB).


## RSCU_for_many_files.py

Call  create_RSCU_heatmap(directiories, path_to_directories) to generate heatmap.
directores - list of names of direcories with files with sequences (pass single directory as one element list)
each directory will be sparete subplot.
path_to_directories = path to direcory containing these direcorie
generates reference set of genes

## get_ref_set.py
generates reference set of genes.

## ramps_tAI.py

generates plots with ramps based on tAI
Function create_plot_of_sum_of_genes() takes list of SeqIO object,
that can be gerated by function create_list_of_records(),
where organism is a file with sequences and len_of_gene is minimal lenght of genes
that you want to study.
Another paramether are:
path, where you pass path of a directory where You want to save generated plot,
name, which is a name for generated file
number_of_ramps, where you define hom many ramps the scripts is to generate
codons, where You are to pass dictionary of codons with their relative_adaptiveness values.
This dictionary can be generatet by passing file with these values to create_list_of_records()

## percent_GC.py

Script to generate ramps with GC precent values.
Call count_gc_for_many_files(path_to_folder, path_dest_folder, name_of_file,  len_of_gene, num_of_ramps, len_of_ramps, step, path_for_overall=None, path_for_ramps=None)
function to generate raps.
path_to_folder is path to folder where files with sequences are located
path_dest_folder is path to folder where you want to save output plot
name_of_file is a name for output file
num of ramps is an integer that says hom many ramps the scripts is to generate
len_of_ramps- how long are the ramps
step- number of nucleotydes that separetes ramps (ex. if lenght of ramp is 3, and step is 3,
ramps are not separated, if lenght of ramp is 3, and step is 4, ramps are separated by 1 nucleotyde)
step must by iteration of 3.

## ramps_cai.py

  generates plots with ramps based on CAI
call  create_sum_of_files_cai function and pass number
of files You want to analyse and path to directory where the files are located

## main.py
"generates ramps with GC percent, RSCU heatmap and ramps with CAI values
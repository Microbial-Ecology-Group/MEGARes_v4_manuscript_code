# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 15:08:33 2025

@author: Jared Young
"""
# version 1.0.1 update notes 
## Fixes annotation heirarchy issue by
## adding accession_ID between group and SNP annotations 
import csv
import os
import argparse

# define main argument 
def main():
    # setting up command-linw interface - updated
    parser = argparse.ArgumentParser(description='Generates a concatenated count matrix and annotation file for the raw files given',
                                     usage='python %(prog)s [-h] [-s save_dir] file_dir fn_struc sample_file')

    parser.add_argument('-s', type=str, default='./', help='Path in which to store the generated files')
    parser.add_argument('file_dir', type=str, help='Path where the raw files are stored')
    parser.add_argument('fn_struc', type=str, help='Naming structure uses in all raw files, e.g., "called_SNPs.best_split_"')
    parser.add_argument('sample_file', type=str, help='Path to the file with the names of all the samples')

    args = parser.parse_args()

    # reading raw count data into a data frame - updated

    data_frame = []

    for fn in os.listdir(args.file_dir):
        if args.fn_struc in fn:
            with open (args.file_dir + fn) as file:
                reader = csv.reader(file, delimiter='\t') # read in SNP counts for parsiing 
                for row in reader: # for row in reader do the following 
                    extras_removed = row[0].split('|')[:5] # fixing RequiresSNPConfirmation
                    ascid = extras_removed[0]
                    asc = extras_removed[0].split('_')[1] # get MEG accession number
                    grp_asc = '-'.join([extras_removed[4]]+[asc]) # merge MEG accession with group with sep by '-'
                    snp_ann = '_'.join([grp_asc]+[row[3],row[2]]) # join group (with MEG accession), reference nucleotide, position to start SNP header
                    joined = '|'.join(extras_removed + [ascid] + [snp_ann]) # join partial SNP header to existing megares annotation, separate with |
                    data_frame.append([joined] + row[5].split(',')) # join new (not complete) snp header with SNP counts - counts sep by "," to find 2nd and 3rd position SNPs
                

    # adjusting data frame structure - updated
    
    for i in range(len(data_frame)):       # for the entire length of the dataframe line
        while len(data_frame[i]) > 2:      # while the line has more than 2 sections (IE: 3-4)
            data_frame.append([data_frame[i][0], data_frame[i][-1]]) # append the dataframe annotation to a new line with last part
            data_frame[i] = data_frame[i][:-1] # remove last part from entry 

    
    for i in range(len(data_frame)): # for row i in the entire data frame 
        data_frame[i] =  [data_frame[i][0]] + data_frame[i][1].split('|') # keep first header column, split second column with SNP info and counts 
        data_frame[i] = [data_frame[i][0], data_frame[i][2]] + data_frame[i][4:] # [2] is alt nucleotide, [4:] is per-sample counts (usually [3:], but 3+1 with header)

        data_frame[i][0] = '_'.join([data_frame[i][0], data_frame[i][1]]) # alt nuc is now position 1, changed from '|' join to '_' join
        data_frame[i] = [data_frame[i][0]] + data_frame[i][2:] # data frame [0] is updated with full SNP annotation, but [1] is still alt nuc, so start [2:] to grab counts, remove alt nuc

    # reading in sample names - updated

    all_samples = []  # make all samples list
    s_add = "S_"      # make s_add variable -> adds S_ in front of sample name if name starts with int
    
    with open(args.sample_file) as file:
        for line in file: # will update samp_names to "file" once I merge formatting with command line code
            line = line.rstrip()
            if "." and "/" in line:
                line = line.split("/")
                line = line[-1].split(".")
                line = line[0]
                if line[0].isdigit():
                    line = s_add + line
                    all_samples.append(line)
                else: 
                    all_samples.append(line)
            elif "." in line:
                line = line.split(".")
                line = line[0]
                if line[0].isdigit():
                    line = s_add + line
                    all_samples.append(line)
                else: 
                    all_samples.append(line)
            elif "/" in line: 
                line = line.split("/")
                line = line[-1]
                if line[0].isdigit():
                    line = s_add + line
                    all_samples.append(line)
                else: 
                    all_samples.append(line)
            else:
                if line[0].isdigit():
                    line = s_add + line
                    all_samples.append(line)
                else: 
                    all_samples.append(line)

    # adding title row to data frame - updated

    data_frame.insert(0, ['SNP_accession'] + all_samples)

    # writing final count matrix - updated
    
    with open(args.s + 'utd_ARG_SNP_count_matrix.csv', 'w', newline='') as file:  # space of blank rows between each row will appear without newline = ''
        writer = csv.writer(file)
        writer.writerows(data_frame)

    # generating annotation file data frame - updated
    
    ann_data_frame = []

    for row in data_frame[1:]:
        ann_data_frame.append([row[0], *row[0].split('|')[1:]]) # keep all levels except leading MEG id 

    title_row = ['SNP_accession', 'Type', 'Class', 'Mechanism', 'Group', "ARG_accession", 'SNP']

    ann_data_frame.insert(0, title_row)
    
    # writing count matrix annotation file - update
    with open(args.s + 'utd_ARG_SNP_annotations.csv','w', newline='') as file: # for testing 
        writer = csv.writer(file)
        writer.writerows(ann_data_frame)
    # done

if __name__ == '__main__':
    main()

# test 
#main(snpCaller/ called_SNPs all_samples)
# python3 clean_metaSNV.py ../../youn2635/four_treatment_SNP/output_dir/snpCaller/ called_SNPs.best_split_ ../../youn2635/four_treatment_SNP/all_samples
## Notes
# when sep by "|"
# Notes: metaSNV input structre example: MEG_4223|Drugs|Class|Mech|Group|RequiresSNPConfirmation	-	28	G	9|24|276|136|...	3385|A|.|0|7|69|25|...
#  ... where:                               [0]                                                     [1] [2] [3] [4]                 [5]
# 0 = megares annotation 
# 1 = either gene containing positioin (for chromosomes) or "-"
# 2 = position within the reference 
# 3 = reference allele (nucleotide)
# 4 = reference allele frequency string (per-read coverage at a specific position in reference), count per sample (in order of input) sep by "|"
# 5 = alternate allele (SNP) frequency string - abundance_all_samples|SNP_nuc|.|x1|x2|x3|x4|...
#                                               [0]                   [1]     [2] [3]
## when line.sep[5] = sep by "|"
# 0 = abundance of SNP in all samples (total reads containing SNP)
# 1 = alternate nucleotide (ie, the SNP)
# 2 = codon change or "-"
# 3 = allele freqency string - total reads per sample containing that SNP, separated by "|", in order of sample input
## Note, more than one alt nucleotide can appear at the same position, separate by a ","
# if line.sep[5] by "," -> abund|A|.|1|2|3|4|...,abund|C|.|1|2|3|4|...,abund|T|.|1|2|3|4|...
#                          [0]                   [1]                   [2]
# 0 = alternate allele string 1 (always there)
# 1 = alternate allele string 2 (may or may not be present)
# 2 = alternate allele string 3 (may or may not be present, least likely to be seen, but does happen)

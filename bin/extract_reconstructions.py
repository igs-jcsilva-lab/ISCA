#!/usr/bin/env python3

'''
This script follows threaded_assess_alignment.py and pulls the longest assembled sequence
regardless of the threshold for each locus into a FASTA file. The FASTA header shall
include the locus ID, %ID, length of the sequence assembled and proportion of reference 
sequence length that aligned with assembled sequence.

This script runs after both the iterations are completed (after first and second_ids_v_cov.tsv generated)
or after the first iteration in cases where second iteration with SMALT isn't necessary.

    Input:
        1. A base directory location
        2. Path to ids_v_cov.tsv file generated from threaded_assess_alignment.py after HGA+SB step in first iteration (first_final_phase_three)
        3. Path to the directory preceding all the alignment directories for first iteration (e.g. for "/path/to/ref123" put "/path/to" as the input)
        4. Path to the directory preceding all the alignment directories for second iteration (e.g. for "/path/to/ref123" put "/path/to" as the input)
        5. Path to ids_v_cov.tsv file generated from threaded_assess_alignment.py after HGA+SB step in second iteration (first_final_phase_three)
        6. Path to ids_v_cov.tsv file generated from threaded_assess_alignment.py after SPAdes step in first iteration (first_intermediary_phase_three)
        7. Path to ids_v_cov.tsv file generated from threaded_assess_alignment.py after SPAdes step in second iteration (second_intermediary_phase_three)
        8. Path to the unbuffered FASTA from extract_sequences.py

    Output:
        An outfile assembled_seqs_below_threshold.fsa containing longest assembled sequences

    Usage:
        extract_reconstructions.py -fivc path/to/file -fap path/to/first_alignments -sap path/to/second_alignments -sivc path/to/file -fidvc path/to/file -sidvc path/to/file

    Author:
        Chakshu Gandhi
'''

import re, os, argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def main():

    parser = argparse.ArgumentParser(description="Script for extracting reconstructions that don't meet threshold")
    parser.add_argument('--workspace_location', '-wl', type=str, default='.', help='Path where directories are built at.')
    parser.add_argument('--first_ids_v_cov', '-fivc', type=str,  help='Path to .tsv file generated from threaded_assess_alignment.py in first iteration after HGA+SB step')
    parser.add_argument('--first_align_path', '-fap', type=str, required=True, help='Path to the directory preceding all the alignment directories for first iteration (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('--second_align_path', '-sap', type=str, required=True, help='Path to the directory preceding all the alignment directories for second iteration (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('--second_ids_v_cov', '-sivc', type=str, help='Path to .tsv file generated from threaded_assess_alignment.py in second iteration after HGA+SB step')
    parser.add_argument('--first_idvc', '-fidvc', type=str,  help='Path to .tsv file generated from threaded_assess_alignment.py in first iteration after SPAdes step')
    parser.add_argument('--second_idvc', '-sidvc', type=str, help='Path to .tsv file generated from threaded_assess_alignment.py in second iteration after SPAdes step')
    parser.add_argument('--original_fsa', '-of', type=str, required=True, help='Path to where the unbuffered FASTA from extract_sequences.py is.')
    args = parser.parse_args()
    
    sequences_list = []
    first_ivc = args.first_ids_v_cov
    second_ivc = args.second_ids_v_cov
    first_idvc = args.first_idvc
    second_idvc = args.second_idvc
    first_ap = args.first_align_path
    second_ap = args.second_align_path
    best_len = defaultdict(list)
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.original_fsa,"fasta"))
    outfile = "{0}/assembled_seqs_below_threshold.fsa".format(args.workspace_location)

    best_len = extract_list(first_idvc, best_len, first_ap, second_ap)
    best_len = extract_list(first_ivc, best_len, first_ap, second_ap)
    best_len = extract_list(second_idvc, best_len, first_ap, second_ap)
    best_len = extract_list(second_ivc, best_len, first_ap, second_ap)
    
    
    # Pulling FASTA sequences together, editing headers for FASTA records
    for k in best_len:
        path = best_len[k][5]
        path_list = path.split('/')
        ref_seq = path_list[-1].split('.WITH')[0]
        new_id = "assembled_{0}".format(k)
        file_ = path.replace('trimmed_align.txt', 'b.fsa')
        
        ref_len = int(best_len[k][4])
        seq = seq_dict[ref_seq]
        ref_len_percent = round((ref_len/len(seq)*100), 2)

        record = SeqIO.read(file_,"fasta")
        record.id = new_id
        record.description = 'ID_to_ref={0} len={1} ref_len_percent={2}'.format(str(best_len[k][0]),str(best_len[k][2]),str(ref_len_percent))
        record.name = ''
        if '.r.trimmed' in path:
            rev_seq = Seq(str(record.seq))
            record.seq = rev_seq.reverse_complement()
        sequences_list.append(record)

    SeqIO.write(sequences_list, outfile, 'fasta')

# Parses through each ids_v_cov.tsv (first and second) and stores the longest assembled
# sequence details in best_len dictionary 
def extract_list(infile, best_len, first_ap, second_ap):
  if os.path.exists(infile):
    with open(infile, 'r') as f1:
        for line in f1:
            line = line.rstrip()
            ivc_list = line.split()
            aligned = ivc_list[5].split('.WITH.')[1]
            locus = ivc_list[5].split('.')[1]
            length = int(ivc_list[4])
            align_dir = ivc_list[5].split('/')[-3]
           
           # Checks if the ids_v_cov.tsv file has paths from first/second_alignments directory
           # and accordingly chooses the align_path to extract sequences
            first_split_point = os.path.basename(first_ap)
            second_split_point = os.path.basename(second_ap)
            if align_dir == first_split_point:
                tmp_path = ivc_list[5].split(first_split_point)[1]
                file_path = "{0}/{1}".format(first_ap,tmp_path)
            elif align_dir == second_split_point:
                tmp_path = ivc_list[5].split(second_split_point)[1]
                file_path = "{0}/{1}".format(second_ap,tmp_path)
            ivc_list[5] = file_path
                
           # Updating best_len dictionary with longest sequence details for each locus     
            if locus in best_len:
                if length > int(best_len[locus][4]) and len(ivc_list) == 6:
                    best_len[locus][0] = float(ivc_list[0])
                    best_len[locus][1] = float(ivc_list[1])
                    best_len[locus][2] = int(ivc_list[2])
                    best_len[locus][3] = ivc_list[3]
                    best_len[locus][4] = int(ivc_list[4])
                    best_len[locus][5] = ivc_list[5]
                elif length > int(best_len[locus][4]) and len(ivc_list) != 6:
                    best_len[locus][0] = float(ivc_list[6])
                    best_len[locus][1] = float(ivc_list[1])
                    best_len[locus][2] = int(ivc_list[2])
                    best_len[locus][3] = ivc_list[3]
                    best_len[locus][4] = int(ivc_list[4])
                    best_len[locus][5] = ivc_list[5]
                elif length == int(best_len[locus][4]) and len(ivc_list) == 6:
                    if float(ivc_list[0]) > best_len[locus][0]:
                        best_len[locus][0] = float(ivc_list[0])
                        best_len[locus][1] = float(ivc_list[1])
                        best_len[locus][2] = int(ivc_list[2])
                        best_len[locus][3] = ivc_list[3]
                        best_len[locus][4] = int(ivc_list[4])
                        best_len[locus][5] = ivc_list[5]
                elif length == int(best_len[locus][4]) and len(ivc_list) != 6:
                    if float(ivc_list[6]) > best_len[locus][0]:
                        best_len[locus][0] = float(ivc_list[0])
                        best_len[locus][1] = float(ivc_list[1])
                        best_len[locus][2] = int(ivc_list[2])
                        best_len[locus][3] = ivc_list[3]
                        best_len[locus][4] = int(ivc_list[4])
                        best_len[locus][5] = ivc_list[5]
                else:
                    continue
            else:
                if len(ivc_list) == 6:
                    best_len[locus].append(float(ivc_list[0]))
                    best_len[locus].append(float(ivc_list[1]))
                    best_len[locus].append(int(ivc_list[2]))
                    best_len[locus].append(ivc_list[3])
                    best_len[locus].append(int(ivc_list[4]))
                    best_len[locus].append(ivc_list[5])
                else:
                    best_len[locus].append(float(ivc_list[6]))
                    best_len[locus].append(float(ivc_list[1]))
                    best_len[locus].append(int(ivc_list[2]))
                    best_len[locus].append(ivc_list[3])
                    best_len[locus].append(int(ivc_list[4]))
                    best_len[locus].append(ivc_list[5])

    return best_len

main()

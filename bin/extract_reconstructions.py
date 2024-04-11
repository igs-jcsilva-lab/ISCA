#!/usr/bin/env python3

'''
This script follows threaded_assess_alignment.py and pulls the longest assembled sequence
regardless of the threshold for each locus and the best %ID sequence for each locus into 
a FASTA file. This script ignores repeated sequences by comparing with assembled_seqs.fsa
file. The FASTA header shall include the locus ID, %ID, length of the sequence assembled and 
proportion of reference sequence length that aligned with assembled sequence.

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
        8. Path to assembled_seqs.fsa file generated from gather_sequences step 
        9. Path to the unbuffered FASTA from extract_sequences.py

    Output:
        1. A file secondary_seqs_by_length.fsa containing longest assembled sequences
        2. Another file secondary_seqs_by_id.fsa containing best %id sequences

    Usage:
        extract_reconstructions.py -fivc path/to/file -fap path/to/first_alignments -sap path/to/second_alignments -sivc path/to/file -fidvc path/to/file -sidvc path/to/file -af path/to/assembled_seqs.fsa -of path/to/unbuffered/fasta/file

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
    parser.add_argument('--assmb_file', '-af', type=str, help='Path to assembled_seqs.fsa file generated from gather_sequences step')
    parser.add_argument('--original_fsa', '-of', type=str, required=True, help='Path to where the unbuffered FASTA from extract_sequences.py is.')
    args = parser.parse_args()
    
    best_len_seq_list = []
    best_id_seq_list = []
    first_ivc = args.first_ids_v_cov
    second_ivc = args.second_ids_v_cov
    first_idvc = args.first_idvc
    second_idvc = args.second_idvc
    first_ap = args.first_align_path
    second_ap = args.second_align_path
    best_len = defaultdict(list)
    best_id = defaultdict(list)
    primary_seqs = defaultdict(list)
    assmb_file = args.assmb_file
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.original_fsa,"fasta"))
    outfile = "{0}/secondary_seqs_by_length.fsa".format(args.workspace_location)
    outfile2 = "{0}/secondary_seqs_by_id.fsa".format(args.workspace_location)
    primary_seqs = primary_output_to_dict(assmb_file,primary_seqs)
    best_len = extract_length_list(first_idvc, best_len, first_ap, second_ap, primary_seqs)
    best_len = extract_length_list(first_ivc, best_len, first_ap, second_ap, primary_seqs)
    best_len = extract_length_list(second_idvc, best_len, first_ap, second_ap, primary_seqs)
    best_len = extract_length_list(second_ivc, best_len, first_ap, second_ap, primary_seqs)
    best_id = extract_id_list(first_idvc, best_id, first_ap, second_ap, primary_seqs)
    best_id = extract_id_list(first_ivc, best_id, first_ap, second_ap, primary_seqs)
    best_id = extract_id_list(second_idvc, best_id, first_ap, second_ap, primary_seqs)
    best_id = extract_id_list(second_ivc, best_id, first_ap, second_ap, primary_seqs)

    
    # Pulling FASTA sequences with best length together, editing headers for FASTA records
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
        best_len_seq_list.append(record)

    SeqIO.write(best_len_seq_list, outfile, 'fasta')

    # Pulling FASTA sequences with best id together, editing headers for FASTA records
    for k in best_id:
        path = best_id[k][5]
        path_list = path.split('/')
        ref_seq = path_list[-1].split('.WITH')[0]
        new_id = "assembled_{0}".format(k)
        file_ = path.replace('trimmed_align.txt', 'b.fsa')

        ref_len = int(best_id[k][4])
        seq = seq_dict[ref_seq]
        ref_len_percent = round((ref_len/len(seq)*100), 2)

        record = SeqIO.read(file_,"fasta")
        record.id = new_id
        record.description = 'ID_to_ref={0} len={1} ref_len_percent={2}'.format(str(best_id[k][0]),str(best_id[k][2]),str(ref_len_percent))
        record.name = ''
        if '.r.trimmed' in path:
            rev_seq = Seq(str(record.seq))
            record.seq = rev_seq.reverse_complement()
        
        best_id_seq_list.append(record)

    SeqIO.write(best_id_seq_list, outfile2, 'fasta')    

# Parses through each ids_v_cov.tsv (first and second) and stores the longest assembled
# sequence details in best_len dictionary 
def extract_length_list(infile, best_len, first_ap, second_ap, primary_seqs):
  if os.path.exists(infile):
    with open(infile, 'r') as f1:
        for line in f1:
            line = line.rstrip()
            ivc_list = line.split()
            aligned = ivc_list[5].split('.WITH.')[1]
            locus = ivc_list[5].split('.')[1]
            length = int(ivc_list[2])
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


            # Ensuring sequences from assembled_seqs.fsa are not repeated
            if primary_seqs[locus]:
                # check length id
                if (length == primary_seqs[locus][1]) and ((len(ivc_list)==6 and float(ivc_list[0])==primary_seqs[locus][0]) or (len(ivc_list) != 6 and float(ivc_list[6])==primary_seqs[locus][0])):
                        continue


           # Updating best_len dictionary with longest sequence details for each locus     
            if locus in best_len:
                if length > int(best_len[locus][2]) and len(ivc_list) == 6 :
                    best_len[locus][0] = float(ivc_list[0])
                    best_len[locus][1] = float(ivc_list[1])
                    best_len[locus][2] = int(ivc_list[2])
                    best_len[locus][3] = ivc_list[3]
                    best_len[locus][4] = int(ivc_list[4])
                    best_len[locus][5] = ivc_list[5]
                elif length > int(best_len[locus][2]) and len(ivc_list) != 6:
                    best_len[locus][0] = float(ivc_list[6])
                    best_len[locus][1] = float(ivc_list[1])
                    best_len[locus][2] = int(ivc_list[2])
                    best_len[locus][3] = ivc_list[3]
                    best_len[locus][4] = int(ivc_list[4])
                    best_len[locus][5] = ivc_list[5]
                elif length == int(best_len[locus][2]) and len(ivc_list) == 6:
                    if float(ivc_list[0]) > best_len[locus][0]:
                        best_len[locus][0] = float(ivc_list[0])
                        best_len[locus][1] = float(ivc_list[1])
                        best_len[locus][2] = int(ivc_list[2])
                        best_len[locus][3] = ivc_list[3]
                        best_len[locus][4] = int(ivc_list[4])
                        best_len[locus][5] = ivc_list[5]
                elif length == int(best_len[locus][2]) and len(ivc_list) != 6:
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

# Parses through each ids_v_cov.tsv (first and second) and stores the best id
# sequence details in best_id dictionary 
def extract_id_list(infile, best_id, first_ap, second_ap, primary_seqs):
  if os.path.exists(infile):
    with open(infile, 'r') as f1:
        for line in f1:
            line = line.rstrip()
            ivc_list = line.split()
            aligned = ivc_list[5].split('.WITH.')[1]
            locus = ivc_list[5].split('.')[1]
            prcnt_id = float(ivc_list[0]) if len(ivc_list)==6 else float(ivc_list[6])
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


            # Ensuring sequences from assembled_seqs.fsa are not repeated
            if primary_seqs[locus]:
                # check length and id
                if (int(ivc_list[2]) == primary_seqs[locus][1]) and ((len(ivc_list)==6 and float(ivc_list[0])==primary_seqs[locus][0]) or (len(ivc_list) != 6 and float(ivc_list[6])==primary_seqs[locus][0])):
                        continue


           # Updating best_id dictionary with best percent id sequence details for each locus     
            if locus in best_id:
                if prcnt_id > float(best_id[locus][0]) and len(ivc_list) == 6 :
                    best_id[locus][0] = float(ivc_list[0])
                    best_id[locus][1] = float(ivc_list[1])
                    best_id[locus][2] = int(ivc_list[2])
                    best_id[locus][3] = ivc_list[3]
                    best_id[locus][4] = int(ivc_list[4])
                    best_id[locus][5] = ivc_list[5]
                elif prcnt_id > float(best_id[locus][0]) and len(ivc_list) != 6:
                    best_id[locus][0] = float(ivc_list[6])
                    best_id[locus][1] = float(ivc_list[1])
                    best_id[locus][2] = int(ivc_list[2])
                    best_id[locus][3] = ivc_list[3]
                    best_id[locus][4] = int(ivc_list[4])
                    best_id[locus][5] = ivc_list[5]
                elif prcnt_id == int(best_id[locus][0]) and len(ivc_list) == 6:
                    if float(ivc_list[0]) > best_id[locus][0]:
                        best_id[locus][0] = float(ivc_list[0])
                        best_id[locus][1] = float(ivc_list[1])
                        best_id[locus][2] = int(ivc_list[2])
                        best_id[locus][3] = ivc_list[3]
                        best_id[locus][4] = int(ivc_list[4])
                        best_id[locus][5] = ivc_list[5]
                elif prcnt_id == int(best_id[locus][0]) and len(ivc_list) != 6:
                    if float(ivc_list[6]) > best_id[locus][0]:
                        best_id[locus][0] = float(ivc_list[0])
                        best_id[locus][1] = float(ivc_list[1])
                        best_id[locus][2] = int(ivc_list[2])
                        best_id[locus][3] = ivc_list[3]
                        best_id[locus][4] = int(ivc_list[4])
                        best_id[locus][5] = ivc_list[5]
                else:
                    continue
            else:
                if len(ivc_list) == 6:
                    best_id[locus].append(float(ivc_list[0]))
                    best_id[locus].append(float(ivc_list[1]))
                    best_id[locus].append(int(ivc_list[2]))
                    best_id[locus].append(ivc_list[3])
                    best_id[locus].append(int(ivc_list[4]))
                    best_id[locus].append(ivc_list[5])
                else:
                    best_id[locus].append(float(ivc_list[6]))
                    best_id[locus].append(float(ivc_list[1]))
                    best_id[locus].append(int(ivc_list[2]))
                    best_id[locus].append(ivc_list[3])
                    best_id[locus].append(int(ivc_list[4]))
                    best_id[locus].append(ivc_list[5])

    return best_id

# Stores the percent id and sequence length details from assembled_seqs.fsa
# in a dictionary for easy comparison
def primary_output_to_dict(assmb_file,primary_seqs):
    record_id_regex = r"[>a-zA-Z]*[_]([a-zA-Z0-9]*[_][0-9]*)"
    prcnt_id_regex = r"[A-Z]+[_][a-z]+[_][a-z]+[=]([0-9.]*)"
    record_len_regex = r"[0-9.]*\s[a-z]+[=]([0-9]+)"
    for record in SeqIO.parse(assmb_file,'fasta'):
        record_id = re.search(record_id_regex,record.id).group(1)
        prcnt_id = re.search(prcnt_id_regex,record.description).group(1)
        record_len = re.search(record_len_regex,record.description).group(1)
        primary_seqs[record_id].append(float(prcnt_id))
        primary_seqs[record_id].append(int(record_len))
    return primary_seqs

main()

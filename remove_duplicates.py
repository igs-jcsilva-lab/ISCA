

# This script goes through a FASTA file and removes subsequent duplicates
# from the file so that only one instance of a given gene sequence is 
# left in the final file. It will also output to STDOUT exactly how many
# duplicates were removed and if any of these sequences were found across
# different loci. Note that the order of the input matters as this script
# will keep the first instance of a duplicate encountered starting from
# the top of the file. This should not be too much of an issue if 
# the input for extract_alleles.py was formatted correctly to guarantee
# that the reference provides the first set of sequences to the file. 
# The number of duplicates output are those that were removed and had
# consistency across loci (so removing dupes found by GMAP). The number
# of conflicts output are those sequences which match but have different
# loci. Take note of this as it may be significant in revealing 
# necessary redundancy in the genome and duplicates here may have been
# removed that should not have. 
#
# *** NOTE ***
# If the number of dupes is non-zero then the follow-up script to this, 
# mod_ea_map.py must be run to purge all duplicate entries from that file.
#
# Run the script using a command like this:
# python3 remove_duplicates.py -input in.fsa -output o.fsa -conflicts conflict.out
#
# Author: James Matsumura

import argparse,re,hashlib
from Bio import SeqIO
from collections import defaultdict
from shared_fxns import write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to refine a FASTA file to not have duplicate sequences.')
    parser.add_argument('-input', type=str, required=True, help='Path to input FASTA.')
    parser.add_argument('-output', type=str, required=True, help='Path to output FASTA.')
    parser.add_argument('-conflicts', type=str, required=True, help='Output of out where conflicts were found.')
    args = parser.parse_args()

    unique_dict = {} # keys are seq hashes and values are locus IDs
    duplicate_entries = []
    duplicates,conflicts = (0 for i in range(2))

    # Just need to iterate through once since we've already preferred the 
    # the reference as the first sequence in extract_sequences.py.
    for record in SeqIO.parse(args.input,"fasta"):
        id = record.id
        locus = re.search(r'\.([A-Za-z0-9_]+)\.?\d?',id).group(1)
        seq = str(record.seq)
        md5_seq = hashlib.md5(seq.encode('utf-8')).hexdigest()
        # Found the first instance of a sequence, write it out
        if md5_seq not in unique_dict:
            unique_dict[md5_seq] = locus
            write_fasta(args.output,id,seq)
        else:
            # If we find an entry with the same locus, remove this new one
            if locus == unique_dict[md5_seq]:
                duplicates += 1
            else: # Same sequence but different loci
                conflicts += 1

            duplicate_entries = process_duplicate_entry(duplicate_entries,id)

    print("Number of duplicates removed: {0}".format(duplicates))
    print("Number of conflicts sequences removed: {0}".format(conflicts))

    with open(args.conflicts,'w') as outfile:
        for dupe in duplicate_entries:
            outfile.write("{0}\n".format(dupe))

# Function to add either just the locus or the full path of the sequence.
# e.g., it converts PF.123456 and PF.123456.1 into 123456 and PF.123456.1 
# respectively.
# Arguments: 
# dupe_list - list to build upon duplicate entries
# id - ID extracted from the header of the FASTA sequence
def process_duplicate_entry(dupe_list,id):

    if re.search(r'\.\d+$',id):
        dupe_list.append(id)
    elif re.search(r'\.path\d+$',id):
        dupe_list.append(id) 
    else:
        dupe_list.append(id.split('.')[1])

    return dupe_list


if __name__ == '__main__':
    main()

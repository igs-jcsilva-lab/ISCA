#!/usr/bin/env python3

"""
A script which follows align_two_sequence_files.py and iterates through all 
directories looking for files ending with the *align.txt in the name. Pulls
the best %ID found for each particular identity and generates a CSV file with 
just two columns, one for the locus name and one for the percent identity. 

Input: 
    1. Path to a directory to iterate through to check alignment stats in. 
    2. Name of the output CSV to generate. 

Output: 
    1. A two-column CSV file with locus_name,percent_id 

Usage:
    assess_all_alignments.py -i /path/to/align/dir -o name_of_out.csv

Author: 
    James Matsumura
"""

import argparse
from collections import defaultdict
import os
import re


def main():

    parser = argparse.ArgumentParser(description='Script to assess alignment results from align_two_sequence_files.py.')
    parser.add_argument('-i', type=str, required=True, help='Path of the base alignment result directory from align_two_sequence_files.py.')
    parser.add_argument('-o', type=str, required=True, help='Name of the CSV file to generate.')
    args = parser.parse_args()

    best_ids = defaultdict(lambda : 0) # map the best ID to each individual locus

    for path,subdirs,files in os.walk(args.i):
        for name in files:
            if name.endswith('align.txt'):

                file_path = "{}/{}".format(path,name)
                locus = file_path.split('/')[-2]

                with open(file_path,'r') as alignment:
                    content = alignment.read()
                    percent_id = float(re.search(r'Identity:\s+\d+/\d+\s\(\s*(\d+\.\d)%\)\n',content).group(1))
                    if percent_id > best_ids[locus]:
                        best_ids[locus] = percent_id

    with open(args.o,'w') as outfile:

        outfile.write("locus_name,percent_id\n")

        for k in sorted(best_ids):
            outfile.write("{},{}\n".format(k,best_ids[k]))


if __name__ == '__main__':
    main()
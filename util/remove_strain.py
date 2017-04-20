

# This script removes a strain from a FASTA file given a prefix that is to
# be ignored. This is built to follow extract_sequences.py.
#
# Run the script using a command like this:
# python3 remove_strain.py -remove 3D7 -original_fsa /path/to/original_fsa -new_fsa /path/to/filtered_fsa
#
# Author: James Matsumura

import re,argparse
from shared_fxns import write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to assess the results of the base Targeted Assembly pipeline.')
    parser.add_argument('-remove', type=str, required=True, help='Prefix of a strain (3D7) to remove or strains, must be CSV: (3D7,7G8).')
    parser.add_argument('-original_fsa', type=str, required=True, help='Path to where the initial FASTA file generated from the pipeline is.')
    parser.add_argument('-new_fsa', type=str, required=True, help='Path to where the output for this filtered FASTA file should go.')
    args = parser.parse_args()

    nonrelevant_alleles = set() # ignore these alleles
    contigs = {} # final FASTA dict

    remove_us = []

    if ',' in args.remove:
        remove_us = args.remove.split(',')
    else:
        remove_us.append(args.remove)

    regex_for_contig_id = ">([a-zA-Z0-9_\.]+)"

    with open(args.original_fsa,'r') as fasta_in:
        for line in fasta_in: # iterate over the FASTA file and extract the entirety of each sequence
            
            line = line.rstrip()

            if line.startswith('>'):

                current_id = re.search(regex_for_contig_id,line).group(1)
                contigs[current_id] = ""

                prefix = line[1:].split('.')[0]

                if prefix in remove_us:
                    nonrelevant_alleles.add(current_id)

            else:
                contigs[current_id] += line # add all the bases

    for allele in contigs:
        if allele not in nonrelevant_alleles:
            write_fasta(args.new_fsa,allele,contigs[allele])

if __name__ == '__main__':
    main()
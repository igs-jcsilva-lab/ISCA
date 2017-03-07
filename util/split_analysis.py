

# This script takes in a FASTA file formatted from the targeted_assembly pipeline
# (primarily after extract_next_round_seqs.py) in order to split the output from
# analyze_sequences.py into two groups: assembled and not assembled.
# a.tsv will be the former and na.tsv will be the latter. 
#
# Run the script using a command like this:
# python3 split_analysis.py -fasta /path/to/inp.fsa -analyze_out /path/to/analyze_sequences.tsv -out_dir output_dir/
#
# Author: James Matsumura

import argparse,os

def main():

    parser = argparse.ArgumentParser(description='Script to assess stats of those sequences that are unaligned.')
    parser.add_argument('-fasta', type=str, required=True, help='Path to FASTA file generated from the targeted assembly pipeline. Must contain reads that cannot assemble.')
    parser.add_argument('-analyze_out', type=str, required=True, help='Path to the TSV file generated from analyze_sequences.py.')
    parser.add_argument('-out_dir', type=str, required=True, help='Path to where the output files should be created.')
    args = parser.parse_args()

    cant_assemble = set()

    with open(args.fasta,'r') as fasta:
        for line in fasta:
            if line.startswith('>'): # only care about headers
                cant_assemble.add(line[1:-1]) # trim ">" and EOL

    # Open up the TSV output from analyze_sequences.py and split these into two groups
    # depending on whether these were assembled or not assembled yet. 
    with open(args.analyze_out,'r') as analysis:
        with open("{0}/a.tsv".format(args.out_dir),'w') as assembled:
            with open("{0}/na.tsv".format(args.out_dir),'w') as not_assembled:
                for line in analysis:
                    id = line.split('\t')[0]

                    if id in cant_assemble:
                        not_assembled.write(line)

                    else:
                        assembled.write(line)

if __name__ == '__main__':
    main()

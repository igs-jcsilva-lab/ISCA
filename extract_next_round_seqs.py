

# This script follows threaded_assess_alignment.py. This is only necessary
# to run if there are still sequences that have yet to be assembled at an 
# adequate threshold. The input is the ids_v_cov.tsv file from the end
# of the pipeline and also the ORIGINAL FASTA file generated. It 
# also requires a minimum threshold for those %ID of the new assembled 
# sequences to have aligned at. This means if you give it a value of 80
# then only %ID alignments from assemblies less than 80 will be put into
# the next round of reference sequences. 
#
# Note that beyond the second round of reconstruction you must concatenate 
# all ids_v_cov.tsv files from previous runs in order to make sure to not
# include a locus that has already been assembled well enough. 
# 
# Run the script using a command like this:
# python3 extract_new_round_seqs.py -ivc /path/to/ids_v_cov.tsv -threshold 80 -original_fsa /path/to/old.fsa -new_fsa /path/to/new.fsa
#
# Author: James Matsumura

import re,argparse
from shared_fxns import write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to map alleles across GFF3 file. Read the top of the file for more details.')
    parser.add_argument('-ivc', type=str, required=True, help='Path to an ids_v_cov.tsv file from the previous run.')
    parser.add_argument('-threshold', type=float, required=True, help='Minimum threshold of %ID that needs to be met to pass final assembly.')
    parser.add_argument('-original_fsa', type=str, required=True, help='Path to where the initial FASTA file generated from the pipeline is.')
    parser.add_argument('-new_fsa', type=str, required=True, help='Path to where the new FASTA file should go.')
    args = parser.parse_args()

    passed = set()
    regex_for_locus = r'/([a-zA-Z0-9\_\.]+).txt'

    # Iterate over the ids_v_cov.tsv file and find those loci which were
    # ABLE to assemble at the minimum threshold. Note that it needs to be
    # done in this manner since not every locus from the original set
    # may have even produced alignments via Bowtie.  
    with open(args.ivc,'r') as i:
        for line in i:
            line = line.rstrip()
            result = line.split('\t') 

            # First check if this sequence passed the minimum threshold
            if float(result[0]) >= args.threshold:

                # Know that the locus is the second group of the alignment.txt
                # file split by periods
                filename = re.search(regex_for_locus,result[3]).group(1)
                locus = filename.split('.')[1]

                passed.add(locus)

    extract_sequences(args.original_fsa,passed,args.new_fsa)


# Arguments:
# file = FASTA file
# loci = set of loci that are already assembled well enough.  
# outfile = output file to write to. 
def extract_sequences(file,loci,outfile):

    regex_for_contig_id = ">([a-zA-Z0-9_\.]+)"

    # Since PF is fairly small, can be greedy about how the FASTA entries are being
    # processed and store them in memory.
    contigs = {}
    keep_us = [] # note which IDs should be re-added at the end
    current_id = "" # store the previous key for the bases to be assigned to

    with open(file,'r') as fasta:
        for line in fasta: # iterate over the FASTA file and extract the entirety of each sequence
            
            line = line.rstrip()

            if line.startswith('>'):

                current_id = re.search(regex_for_contig_id,line).group(1)
                contigs[current_id] = ""

                # in addition to grabbing entire header, check if this entry is needed later
                locus = line.split('.')[1] 
                if locus not in loci:
                    keep_us.append(current_id)

            else:
                contigs[current_id] += line # add all the bases

    for allele in keep_us:

        # Print out in standard FASTA format
        write_fasta(outfile,allele,contigs[allele])


if __name__ == '__main__':
    main()
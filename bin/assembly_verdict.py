

# This script follows threaded_assess_alignment.py. This is only necessary
# to run if there are still sequences that have yet to be assembled at an 
# adequate threshold. The input is the ids_v_cov.tsv file from the end
# of the pipeline, the ORIGINAL FASTA file generated, and the output map 
# from format_for_assembly.py. It also requires a minimum threshold for 
# those %ID of the new assembled sequences to have aligned at. This means 
# if you give it a value of 80 then only %ID alignments from assemblies 
# less than 80 will be put into the next round of reference sequences. 
#
# Note that beyond the second round of reconstruction you must concatenate 
# all ids_v_cov.tsv files from previous runs in order to make sure to not
# include a locus that has already been assembled well enough. 
# 
# The output will be three files. Sequence sets for those that could not align
# and those that could not be assembled in two separate files as well as 
# a new map for assembling. This map can be used by run_HGA.py to run
# a new set of assemblies for these particular loci.  
#
# Run the script using a command like this:
# python3 final_verdict.py --ivc /path/to/ids_v_cov.tsv --threshold 80 --original_assmb_map /path/to/assmb_map.tsv --original_fsa /path/to/old.fsa --out_dir /path/to/new.fsa
#
# Author: James Matsumura

import re,argparse
from shared_fxns import write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to assess the results of the base Targeted Assembly pipeline.')
    parser.add_argument('--ivc', '-i', type=str, required=True, help='Path to an ids_v_cov.tsv file from the previous run.')
    parser.add_argument('--threshold', '-t', type=float, required=True, help='Minimum threshold of %ID that needs to be met to pass final assembly.')
    parser.add_argument('--original_fsa', '-of', type=str, required=True, help='Path to where the initial FASTA file generated from the pipeline is.')
    parser.add_argument('--original_assmb_map', '-oam', type=str, required=True, help='Path to where the output from format_for_assembly.py is located.')
    parser.add_argument('--out_dir', '-o', type=str, required=True, help='Path to where the unaligned/unassembled FASTA entries and the new alignments map should go.')
    args = parser.parse_args()

    assembled,aligned = (set() for i in range(2))
    regex_for_locus = r'/([a-zA-Z0-9\_\.]+).txt'

    # Iterate over the ids_v_cov.tsv file and find those loci which were
    # ABLE to assemble at the minimum threshold. Note that it needs to be
    # done in this manner since not every locus from the original set
    # may have even produced alignments from the recruiter.  
    with open(args.ivc,'r') as i:
        for line in i:
            line = line.rstrip()
            result = line.split('\t') 

            filename = re.search(regex_for_locus,result[3]).group(1)
            locus = filename.split('.')[1]

            percent_id = result[0]

            # If we are using the assessment results from scaffold builder, 
            # get the reference % ID to exclude gaps
            if len(result) == 5:
                percent_id = result[4]

            # First check if this sequence passed the minimum threshold
            if float(percent_id) >= args.threshold:

                # Know that the locus is the second group of the alignment.txt
                # file split by periods
                assembled.add(locus)
                aligned.add(locus)
            
            else:
                aligned.add(locus)

    not_assembled = "{0}/not_assembled.fasta".format(args.out_dir)
    not_aligned = "{0}/not_aligned.fasta".format(args.out_dir)
    leftovers = "{0}/total_leftovers.fasta".format(args.out_dir)
    extract_sequences(args.original_fsa,assembled,not_assembled,aligned,not_aligned,leftovers)

    new_assmb_map = "{0}/new_assmb_map.tsv".format(args.out_dir)
    new_id = 1 # start a counter for new SGE ID for this assembly
    with open(new_assmb_map,'w') as o:
        with open(args.original_assmb_map,'r') as i:
            for line in i:
                line = line.rstrip()
                elements = line.split('\t')
                locus = elements[0]
                sge_id = elements[1]
                if locus not in assembled:
                    o.write("{0}\t{1}\t{2}\n".format(locus,sge_id,new_id))
                    new_id += 1


# Arguments:
# file = FASTA file
# loci1 = set of loci that are already assembled well enough.  
# loci2 = set of loci that at least aligned and tried to assemble.
# outfile1 = output file to write unassembled to. 
# outfile2 = output file to write unaligned to. 
# outfile3 = both contents of outfile1 and outfile2
def extract_sequences(file,assembled,outfile1,aligned,outfile2,outfile3):

    regex_for_contig_id = ">([a-zA-Z0-9_\.]+)"

    # Since PF is fairly small, can be greedy about how the FASTA entries are being
    # processed and store them in memory.
    contigs = {}
    not_assembled,not_aligned = ([] for i in range(2)) # note which IDs should be re-added at the end
    current_id = "" # store the previous key for the bases to be assigned to

    with open(file,'r') as fasta:
        for line in fasta: # iterate over the FASTA file and extract the entirety of each sequence
            
            line = line.rstrip()

            if line.startswith('>'):

                current_id = re.search(regex_for_contig_id,line).group(1)
                contigs[current_id] = ""

                # in addition to grabbing entire header, check if this entry is needed later
                locus = line.split('.')[1] 
                if locus not in aligned: # we know that if it didn't align, couldn't have assembled
                    not_aligned.append(current_id)
                elif locus not in assembled:
                    not_assembled.append(current_id)

            else:
                contigs[current_id] += line # add all the bases

    for allele in not_assembled:
        write_fasta(outfile1,allele,contigs[allele])
        write_fasta(outfile3,allele,contigs[allele])

    for allele in not_aligned:
        write_fasta(outfile2,allele,contigs[allele])
        write_fasta(outfile3,allele,contigs[allele])


if __name__ == '__main__':
    main()
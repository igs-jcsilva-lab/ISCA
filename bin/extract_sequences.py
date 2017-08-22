#!/usr/bin/env python3

# This script follows extract_alleles.py. This expects the output from that 
# script in addition to the same list used as input for that script. The 
# output of this script will be a FASTA file. The final argument of this script,
# -buffer, is optional. This parameter sets a "buffer" region to extend the gene
# positions noted in the GFF3 file a bit. The output of this is meant to build
# a reference tob e used for alignment.
#
# It is recommended that you set a value for the buffer length that makes sense
# depending on your read length (e.g. use 125 for 250bp reads). If you do set
# a buffer, you should run this script twice as you'll want a set of sequences
# that lack a buffer when it comes to the alignment steps.  
# 
# Run the script using a command like this:
# extract_sequences.py --ea_input /path/to/list_input.tsv --ea_map /path/to/out_from_extract_alleles.tsv --prefix /path/to/my_seqs --buffer 20
#
# Author: James Matsumura

import re,argparse
from Bio import SeqIO
from shared_fxns import write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to extract sequences from each allele across GFF3/FASTA files. Read the top of the file for more details.')
    parser.add_argument('--ea_input', '-eai', type=str, required=True, help='Path to a TSV list for references and isolates.')
    parser.add_argument('--ea_map', '-eam', type=str, required=True, help='Path to the output from extract_alleles.py.')
    parser.add_argument('--subset_list', '-sl', type=str, required=True, help='Path to a list of loci to subset by.')
    parser.add_argument('--buffer', '-b', type=int, default=0, required=False, help='How much of a buffer to add to each end of the gene. Defaults to 0.')
    parser.add_argument('--prefix', '-p', type=str, required=True, help='Prefix of the output FASTA file to generate in current or existing directory.')
    args = parser.parse_args()

    extract_us = {}
    subset_by = set()

    with open(args.subset_list,'r') as i:
        for locus in i:
            subset_by.add(locus.rstrip().split('\t')[0])

    # Iterate over the output from extract_alleles.py and build a dict of lists
    # for all the regions from each FASTA file that need to be extracted. 
    with open(args.ea_map,'r') as i:
        for allele in i: 
            allele = allele.rstrip()
            indv_allele = allele.split('\t')

            if args.subset_list:
                locus = indv_allele[0]
                if locus not in subset_by:
                    continue

            for j in range(1,len(indv_allele)): # iterate over all alleles identified
                vals = indv_allele[j].split('|')
                name = vals[4].split('.')[0] 

                if name not in extract_us: # assign this entry to the proper list via key
                    extract_us[name] = []

                extract_us[name].append(indv_allele[j])

    # Iterate over the input list and extract sequences per FASTA file.
    with open(args.ea_input,'r') as l:
        for entry in l:
            entry = entry.rstrip()
            vals = entry.split('\t')
            fasta_file = vals[2]
            gene_list = extract_us[vals[3]]

            extract_sequences(fasta_file,gene_list,args.buffer,args.prefix)


# Arguments:
# file = FASTA file
# genes = list of genes associated with the particular FASTA file. Noted by the
# name in the fourth column of the input list file.
# buffer = amount of bases to pad the surrounding regions of the gene by.  
# prefix = base prefix to write unbuffered/buffered files to 
def extract_sequences(file,genes,buffer,prefix):

    # Since PF is fairly small, can be greedy about how the FASTA entries are being
    # processed and store them in memory.
    contigs = SeqIO.to_dict(SeqIO.parse(file, "fasta"))

    for gene in genes:
        vals = gene.split('|')
        source = vals[0]
        id = vals[4]
        strand = vals[3]
        start = (int(vals[1]) - buffer - 1) # correct for 0-base indexing
        stop = (int(vals[2]) + buffer) 

        # know that start can't be less than 1 due to buffer alteration, 
        # need to handle stop similarly using the sequence length.
        if start < 1:
            start = 0
        if stop > len(contigs[source].seq):
            stop = len(contigs[source].seq)

        # Generate both a buffered and unbuffered sequence. These files will
        # be exactly the same if buffer==0
        sequence1 = str(contigs[source].seq[start:stop].upper())
        sequence2 = str(contigs[source].seq[(int(vals[1]) - 1):int(vals[2])].upper())

        if strand == "-": # if reverse strand, swap the bases
            sequence1 = rev_comp(sequence1)
            sequence2 = rev_comp(sequence2)

        # Print out in standard FASTA format
        write_fasta("{0}_buffered.fsa".format(prefix),id,sequence1)
        write_fasta("{0}_unbuffered.fsa".format(prefix),id,sequence2)

# Function to build the reverse complement of a sequence
def rev_comp(sequence):

    complement = {"A":"T","T":"A","G":"C","C":"G"}
    mod_seq = ""
    for base in sequence:
        if base in complement:
            mod_seq += complement[base]
        else: # just add Ns
            mod_seq += base

    sequence = mod_seq
    sequence = sequence[::-1]

    return sequence


if __name__ == '__main__':
    main()
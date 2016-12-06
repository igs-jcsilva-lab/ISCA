

# This script follows extract_alleles.py. This expects the output from that 
# script in addition to the same list used as input for that script. The 
# output of this script will be a FASTA file. The final argument of this script,
# -b, is optional. This parameter sets a "buffer" region to extend the gene
# positions noted in the GFF3 file a bit. 
# 
# Run the script using a command like this:
# python3 extract_sequences.py -l /path/to/list_input.tsv -i /path/to/out_from_extract_alleles.tsv -o /path/to/out.fsa -b 20
#
# Author: James Matsumura

import sys,re,argparse

def main():

    parser = argparse.ArgumentParser(description='Script to map alleles across GFF3 file. Read the top of the file for more details.')
    parser.add_argument('-l', type=str, required=True, help='Path to a TSV list for references and isolates.')
    parser.add_argument('-i', type=str, required=True, help='Path to the output from extract_alleles.py.')
    parser.add_argument('-b', type=str, required=False, help='How much of a buffer to add to each end of the gene. Defaults to 0.')
    parser.add_argument('-o', type=str, required=True, help='Path to where the output TSV should go.')
    args = parser.parse_args()

    l = open(args.l,'r')
    i = open(args.i,'r')
    o = open(args.o,'w')
    b = 0 # buffer region defaults to 0
    if args.b:
        b = args.b

    extract_us = {}

    # Iterate over the output from extract_alleles.py and build a dict of lists
    # for all the regions from each FASTA file that need to be extracted. 
    for allele in i: 
        allele = allele.rstrip()
        indv_allele = allele.split('\t')

        for j in range(1,len(indv_allele)): # iterate over all alleles identified
            vals = indv_allele[j].split('|')
            name = vals[4].split('.')[0] 

            if name not in extract_us: # assign this entry to the proper list via key
                extract_us[name] = []

            extract_us[name].append(indv_allele[j])

    # Iterate over the input list and extract sequences per FASTA file.
    for entry in l:
        entry = entry.rstrip()
        vals = entry.split('\t')
        fasta_file = vals[2]
        gene_list = extract_us[vals[3]]

        extract_sequences(fasta_file,gene_list,b,o)

# Arguments:
# file = FASTA file
# genes = list of genes associated with the particular FASTA file. Noted by the
# name in the fourth column of the input list file.
# buffer = amount of bases to pad the surrounding regions of the gene by.  
# outfile = output file to write to. 
def extract_sequences(file,genes,buffer,outfile):

    fasta = open(file,'r') 

    regex_for_contig_id = ">([a-zA-Z0-9_]+)"

    # Since PF is fairly small, can be greedy about how the FASTA entries are being
    # processed and store them in memory.
    contigs = {}
    current_id = "" # store the previous key for the bases to be assigned to

    for line in fasta: # iterate over the FASTA file and extract the entirety of each sequence
        
        line = line.rstrip()

        if line.startswith('>'):
            current_id = re.search(regex_for_contig_id,line).group(1)
            contigs[current_id] = ""
        else:
            contigs[current_id] += line # add all the bases

    # Calculate the lengths of each contig just once to make sure that the
    # buffer does not exceed the max length.
    lengths = {} 
    for key in contigs:
        lengths[key] = len(contigs[key])

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
        if stop > lengths[source]:
            stop = lengths[source]

        sequence = contigs[source][start:stop]
        sequence = sequence.upper() # remove any lowercase

        if strand == "-": # if reverse strand, swap the bases
            sequence = sequence[::-1] # reverse
            complement = {"A":"T","T":"A","G":"C","C":"G"}
            mod_seq = ""
            for base in sequence:
                if base in complement:
                    mod_seq += complement[base]
                else: # just add Ns
                    mod_seq += base

        # Print out in standard FASTA format
        outfile.write(">{0}\n".format(id))
        for i in range(0, len(sequence), 60):
            outfile.write(sequence[i:i+60] + "\n")

if __name__ == '__main__':
    main()
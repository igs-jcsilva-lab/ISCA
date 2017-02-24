

# This script takes in a FASTA file formatted from the targeted_assembly pipeline
# and extracts from the FASTA headers exon counts from a particular GFF3 file.
# The output will be a TSV with 5 columns:
#
# 1. ID
# 2. GC content
# 3. Length
# 4. Number of repeats
# 5. Number of exons. 
#
# The program Red (REpeat Detector) is used to generate an *.rpt file from the given
# FASTA sequence to get an idea of how many repeats are present per sequence. 
# http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0654-5#MOESM4
#
# Run the script using a command like this:
# python3 analyze_genes.py -fasta /path/to/inp.fsa -gff3 pf.gff3 -priority 3D7 -red_out out.rpt -outfile output_dir/out.tsv
#
# Author: James Matsumura

import argparse,os,re
from collections import defaultdict
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(description='Script to assess stats of those sequences that are unaligned.')
    parser.add_argument('-fasta', type=str, required=True, help='Path to FASTA file generated from the targeted assembly pipeline.')
    parser.add_argument('-gff3', type=str, required=True, help='Path to the GFF3 annotation to count exons from.')
    parser.add_argument('-priority', type=str, required=True, help='Prefix of the allele to generate stats for, must match prefix in headers of the input FASTA.')
    parser.add_argument('-red_out', type=str, required=True, help='Path to the .rpt file output from Red.')
    parser.add_argument('-outfile', type=str, required=True, help='Path to the output file.')
    args = parser.parse_args()

    loci = [] # list for outputting locus objects once we have all the info needed
    exon_cnts,repeat_cnts = (defaultdict(lambda:0) for i in range(2)) # count number of exons and repeats

    # Extract from the GFF3 file the number of exons per gene.  
    with open(args.gff3,'r') as gff3:

        exon_cnt = 0
        id = ""

        regex_for_id = r'ID=([a-zA-Z0-9_\-\.]+)'
        regex_for_desc = r'description=(.*)'

        for line in gff3:

            if line.startswith('##FASTA'):
                break
            elif line.startswith('#'):
                continue
            else:
                line = line.rstrip()
                ele = line.split('\t')
                if ele[2] == 'gene': 
                    # Now within a gene, reset exon count until next gene
                    exon_cnt = 0
                    attr = ele[8].split(';')
                    for tag in attr:
                        if tag.startswith('ID'):
                            id = "{0}.{1}".format(args.priority,re.search(regex_for_id,tag).group(1))
                if ele[2] == 'exon':
                    exon_cnts[id] += 1      

    # Extract from the Red output *.rpt file the number of repeats
    with open(args.red_out,'r') as rpt:
        
        for line in rpt:
            id = line.split(':')[0]
            repeat_cnts[id[1:]] += 1

    # Extract the sequences from the FASTA file and *STORE IN MEMORY*
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.fasta,"fasta"))

    # Write out a TSV file with the ID, GC content, length, number of
    # repeats, and number of exons. 
    with open(args.outfile,'w') as out:    
        for key in seq_dict:
            # Find all the IDs we want to work with based on priority.
            if key.startswith(args.priority):
                id = key
                gc = calc_gc_content(seq_dict[id].seq)
                length = len(seq_dict[id].seq)
                repeats = repeat_cnts[id]
                exons = exon_cnts[id]
                out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(id,gc,length,repeats,exons))

# Function to calculate GC percentage given a sequence.
# Argument:
# seq = nucleotide sequence to be analyzed
def calc_gc_content(seq):
    bases = {"A":0,"T":0,"C":0,"G":0}
    for base in seq:
        bases[base] += 1

    gc_cont = (100*(bases["G"]+bases["C"])/(bases["A"]+bases["T"]+bases["G"]+bases["C"]))
    return float("{0:.2f}".format(gc_cont))

if __name__ == '__main__':
    main()

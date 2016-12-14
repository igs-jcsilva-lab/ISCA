

# This script uses EMBOSS's needle global alignment tool to perform alignments
# between the assembled output of format_for_spades.py+SPAdes and the initial 
# FASTA set used to build the Bowtie2 reference.  
#
# Run the script using a command like this:
# python3 global_alignment.py -i /path/to/spades_map.tsv -f /path/to/ref.fasta -out /path/to/out_alignment.info
#
# Author: James Matsumura

import sys,re,argparse
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

def main():

    parser = argparse.ArgumentParser(description='Script to generate stats given output from analyze_bam.py and filter a set of paired-end FASTQ reads.')
    parser.add_argument('-refs', type=str, required=True, help='Path to *_ref_map.tsv output from analyze_bam.py.')
    parser.add_argument('-path', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-out', type=str, required=True, help='Path to output map (maps the ref to the SGE ID)).')
    args = parser.parse_args()

    needle = NeedleCommandline(r"/usr/local/packages/EMBOSS-6.3.1/bin/needle",
                                asequence="/local/scratch/matsu_pf_work/pf_example.fasta",
                                bsequence="/local/scratch/matsu_pf_work/spades_out/PF3D7_1478300/contigs.fasta",
                                gapopen=10, gapextend=0.5,outfile=args.out)

    filename = args.out
    format = "emboss"

    for alignment in AlignIO.parse(filename,format):
        for sequence in alignment:
            print(sequence.seq)


if __name__ == '__main__':
    main()

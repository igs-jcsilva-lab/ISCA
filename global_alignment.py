

# This script uses EMBOSS's needle global alignment tool to perform alignments
# between the assembled output of format_for_spades.py+SPAdes and the initial 
# FASTA set used to build the Bowtie2 reference.  
#
# Run the script using a command like this:
# python3 global_alignment.py -i /path/to/spades_map.tsv -f /path/to/ref.fasta -out /path/to/refined_alignment_out.txt
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

    needle_exe = r"/usr/local/packages/EMBOSS-6.3.1/bin/needle"
    first_align = "{0}/{1}".format(args.out,'first_align.txt')
    format = "emboss"

    needle = NeedleCommandline(r"/usr/local/packages/EMBOSS-6.3.1/bin/needle",
                                asequence="/local/scratch/matsu_pf_work/pf_example.fasta",
                                bsequence="/local/scratch/matsu_pf_work/spades_out/PF3D7_1478300/contigs.fasta",
                                gapopen=10,gapextend=0.5,outfile=first_align)

    stdout,stderr = needle()

    a,b = (None for i in range(2))

    for alignment in AlignIO.parse(first_align,format):
        for sequence in alignment:

            if a == None: # grab both sequences, first being the reference seq
                a = sequence.seq
            else: # now grab the assembled seq
                b = sequence.seq            

            # Once two sequences are extracted, refine and align trimming the 
            # outside extended blank sequence.  
            if a != None and b != None:
                refined_align = "{0}/{1}".format(args.out,'trimmed_align.txt')
                a_fsa = "{0}/{1}".format(args.out,'a.fsa')
                b_fsa = "{0}/{1}".format(args.out,'b.fsa')

                seqs = trim_extensions(a,b)
                write_fasta(a_fsa,'a.trimmed',seqs['a'])
                write_fasta(b_fsa,'b.trimmed',seqs['b'])

                needle = NeedleCommandline(r"/usr/local/packages/EMBOSS-6.3.1/bin/needle",
                                asequence=a_fsa,
                                bsequence=b_fsa,
                                gapopen=10, gapextend=0.5,outfile=refined_align)
                stdout,stderr = needle()


# Function to trim the extended blank bases identified from a Needle alignment.
# Note that this trimming just removes the blanks present in the extension on
# the perimeters of the sequence. 
# Arguments:
# a = first sequence (reference extracted FASTA)
# b = second sequence (SPAdes assembled sequence)
def trim_extensions(a,b):

    a = str(a)
    b = str(b)
    curr_length = len(a) # length before trimming left
    a = a.lstrip('-')
    l_trim = curr_length - len(a)

    curr_length = len(a) # length before trimming right    
    a = a.rstrip('-')
    r_trim = (curr_length - len(a)) * -1

    b = b[l_trim:r_trim]
    a = a.replace('-','') # remove embedded gaps, let Needle re-add
    return {'a':a , 'b':b}

# Function to generate new FASTA files after trimming extensions.
# Arguments:
# file: name/path of file to be written
# header: header ID
# seq: sequence to be written
def write_fasta(file,header,seq):
    with open(file,'w') as out:
        out.write(">{0}\n".format(header))
        for j in range(0, len(seq), 60):
            out.write(seq[j:j+60] + "\n")

if __name__ == '__main__':
    main()

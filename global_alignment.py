

# This script uses EMBOSS's needle global alignment tool to perform alignments
# between the assembled output of format_for_spades.py+SPAdes and the initial 
# FASTA set used to build the Bowtie2 reference.  
#
# Run the script using a command like this:
# python3 global_alignment.py -ea_map /path/to/extract_alleles_map.tsv -ffs_map /path/to/format_for_spades.tsv -ref_genome /path/to/ref_genome.fsa -assmb_path -/path/to/spades_out -out /path/to/alignment_out
#
# Author: James Matsumura

import argparse,os,sys
from collections import defaultdict
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from Bio import SeqIO
from shared_fxns import make_directory,write_fasta

def main():

    parser = argparse.ArgumentParser(description='Script to generate EMBOSS Needle alignments given output from format_for_spades.py.')
    parser.add_argument('-ea_map', type=str, required=True, help='Path to map.tsv output from extract_alleles.py.')
    parser.add_argument('-ffs_map', type=str, required=True, help='Path to map.tsv output from format_for_spades.py.')
    parser.add_argument('-ref_genome', type=str, required=True, help='Path to the reference genome file used to build Bowtie2 index.')
    parser.add_argument('-assmb_path', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-out', type=str, required=True, help='Path to output directory for all these alignments.')
    args = parser.parse_args()

    # First, extract the sequences from the reference file and 
    # *STORE IN MEMORY* (careful how big the reference genome used is.
    # We need this to generate small FASTA files for Needle alignment. 
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.ref_genome,"fasta"))

    # In order to access these seqs efficiently, rebuild the DS created 
    # in extract_alleles.py. This is a dictionary where the key is the 
    # shared locus and the value is a list of all the mapped alleles. 
    ref_dict = defaultdict(list)
    with open(args.ea_map,'r') as loc_map:

        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[0]

            for j in range(1,len(ele)):
                allele_info = ele[j].split('|')
                allele = allele_info[4]
                ref_dict[locus].append(allele)

    # Now that we can easily extract the sequences for alignment, iterate over
    # the directory name map file and perform alignments. 
    with open(args.ffs_map,'r') as dir_map:

        for line in dir_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[0] # reference/locus that maps to directory number
            loc_dir = ele[1] # the directory number from SPAdes for grid submission
            out_dir = "{0}/{1}".format(args.out,locus) # alignment output goes here
            make_directory(out_dir)

            # Split out the contigs if more than one is present and have to do
            # alignment of all refs to all contigs.
            contigs = "{0}/{1}/contigs.fasta".format(args.assmb_path,loc_dir)

            # SPAdes is not able to assemble all the reads, often this seems 
            # to be due to low coverage. Output this to STDOUT. 
            if not os.path.isfile(contigs):
                print("{0}\tcould not assemble.".format(locus))
                continue

            # Iterate over each contig assembled by SPAdes
            for record in SeqIO.parse(contigs,"fasta"):

                type_map = {}
                bseq_file = "{0}/{1}.fsa".format(out_dir,record.id)

                # Make individual FASTA files for each contig
                with open(bseq_file,'w') as bfsa:
                    SeqIO.write(record,bfsa,"fasta")

                # Iterate over each distinct ref sequence (or allele) associated
                # with this particular locus. 
                for ref_seq in ref_dict[locus]:
                    seq = seq_dict[ref_seq]
                    aseq_file = "{0}/{1}.fsa".format(out_dir,ref_seq)

                    # Make an individual FASTA file for each allele
                    if not os.path.isfile(aseq_file): # skip if made for previous contig
                        with open(aseq_file,'w') as afsa:
                            SeqIO.write(seq,afsa,"fasta")

                    # Now have the reference FASTA file, perform alignment
                    # with the assembled contig.
                    alignment = align(out_dir,ref_seq,record.id,aseq_file,bseq_file)
                    align_id = "{0}/{1}.WITH.{2}.trimmed_align.txt".format(out_dir,ref_seq,record.id)
                    type_map[align_id] = alignment['type']

                # For each contig, note the type of alignment so that when the
                # best score is extracted we know whether this could be 
                # re-aligned and potentially better optimized. 
                type_file = "{0}/ga_map.tsv".format(args.out)
                with open(type_file,'a') as type_out:
                    for k,v in type_map.items():
                        type_out.write("{0}\t{1}\n".format(k,v))


# Function to perform an alignment between the reference FASTA sequence and
# the SPAdes contig assembled. 
# Arguments:
# out = output directory for this script
# allele = name of the allele used in the alignment
# aseq = newly generated FASTA file extracted from this script
# bseq = location of the contig.fasta file generated by SPAdes
def align(out,allele,contig,aseq,bseq):

    needle_exe = r"/usr/local/packages/EMBOSS-6.3.1/bin/needle" # path to EMBOSS Needle executable
    format = "emboss"

    initial_align = "{0}/{1}.WITH.{2}.align.txt".format(out,allele,contig)

    needle = NeedleCommandline(needle_exe,
                                asequence=aseq,
                                bsequence=bseq,
                                gapopen=10,gapextend=0.5,outfile=initial_align)
    stdout,stderr = needle()

    a,b = (None for i in range(2))

    # This should just be one alignment... but handles more. 
    for alignment in AlignIO.parse(initial_align,format):
        for sequence in alignment:

            if a == None: # grab both sequences, first being the reference seq
                a = sequence.seq
            else: # now grab the assembled seq
                b = sequence.seq            

            # Once two sequences are extracted, refine and align trimming the 
            # outside extended blank sequence.  
            if a != None and b != None:
                refined_align = "{0}/{1}.WITH.{2}.trimmed_align.txt".format(out,allele,contig)
                a_fsa = "{0}/{1}.WITH.{2}.a.fsa".format(out,allele,contig) # filename
                b_fsa = "{0}/{1}.WITH.{2}.b.fsa".format(out,allele,contig)
                a_trim = "a.trimmed".format(allele,contig) # sequence header, file name makes distinction
                b_trim = "b.trimmed".format(allele,contig)

                seqs = trim_extensions(a,b)
                write_fasta(a_fsa,a_trim,seqs['a'])
                write_fasta(b_fsa,b_trim,seqs['b'])

                needle = NeedleCommandline(needle_exe,
                                asequence=a_fsa,
                                bsequence=b_fsa,
                                gapopen=10, gapextend=0.5,outfile=refined_align)
                stdout,stderr = needle()

                return seqs

# Function to trim the extended blank bases identified from a Needle alignment.
# Note that this trimming just removes the blanks present in the extension on
# the perimeters of the reference locus sequence. If the contig generated by
# SPAdes is shorter than the reference, then no adjustment is made and the
# original alignment should be the same as the trimmed alignment. This returns
# a dictionary of the two sequences as well as the type of alignment found.
# Arguments:
# a = first sequence (reference extracted FASTA)
# b = second sequence (SPAdes assembled sequence)
def trim_extensions(a,b):

    a = str(a)
    b = str(b)

    # check the length of the actual sequence, no gaps!
    aseq = a.replace('-','')
    bseq = b.replace('-','')

    # Leave early if the assembled sequence is smaller than the reference
    # and do not trim.                      ref:      ========
    #                                       assembled:   ===
    if len(aseq) > len(bseq):
        return {'a':aseq,'b':bseq,'type':'ref'}

    # If the sequences overlap like this:   ref:      =========
    #                                       assembled:      ==========
    # it's a bit harder to trim. Seems inappropriate to just chop off the 
    # regions that overextend beyond the alignment. Do not modify these
    # cases for now. 
    if (a[0] == "-" and b[0] != "-") or (b[0] == "-" and a[0] != "-"):
        return {'a':aseq,'b':bseq,'type':'staggered'}

    # If we've made it here, know that the reference likely falls 
    # entirely within the assembly. Trim the overextensions from 
    # the assembled contig. Another case that his point is when 
    # the assembled sequence covers a great deal of the internal 
    # region. Thus, while there is no outside to trim, the assembly
    # is indeed longer than the reference. While unlikely, it may 
    # also be the case where the two are the exact same length. 
    #                                       ref:         ===
    #                                       assembled: =======
    curr_length = len(a) # length before trimming left
    a = a.lstrip('-')
    l_trim = curr_length - len(a)

    curr_length = len(a) # length before trimming right    
    a = a.rstrip('-')
    r_trim = (curr_length - len(a)) * -1

    # If it's the case where the assembly fills in the internal regions,
    # make sure not to do an subset of [0:0] for b. 
    if l_trim == 0 and r_trim == 0:
        pass
    else:
        b = b[l_trim:r_trim]
    a = a.replace('-','') # remove embedded gaps, let Needle re-add
    return {'a':a ,'b':b,'type':'assmb'}


if __name__ == '__main__':
    main()

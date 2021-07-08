#!/usr/bin/env python3

"""
A script which aligns like sequences in two FASTA files. This is meant to help
assess assembled sequences constructed or mapped outside of this pipeline. It
will parse through FASTA headers, split IDs by periods, and align any that match
across files. For instance, ABC.123.path1 and DEF.123.path2 will align to one 
another since '123' match. Each distinct alignment will get its own file in 
the specified output directory. It will load these FASTA files into memory so
careful with their size. 

Input: 
    1. Paths to two FASTA files separated by commas. The first should be a 
    reference file while the second is a set of assembled sequences. 
    2. Path to install directory of EMBOSS needle executable
    3. Location to write all alignments to 

Output: 
    1. Numerous alignment files located in the specified output directory which
    each have names generated from combining the two FASTA header IDs of the 
    two aligned sequences. 

Usage:
    align_two_sequence_files.py -f fasta1.fsa,fasta2.fsa -n /path/to/emboss/bin/needle -o /path/to/out/dir

Author: 
    James Matsumura
"""

import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Emboss.Applications import NeedleCommandline
from shared_fxns import make_directory,write_fasta


def main():

    parser = argparse.ArgumentParser(description='Script to perform needle alignment on like sequences across FASTA files.')
    parser.add_argument('-f', type=str, required=True, help='Two paths to FASTA files split by a comma.')
    parser.add_argument('-n', type=str, required=True, help='Path to install directory of EMBOSS needle executable (e.g. /path/to/packages/emboss/bin/needle).')
    parser.add_argument('-o', type=str, required=True, help='Location to generate output directories.')
    args = parser.parse_args()

    make_directory(args.o)

    first_seqs = SeqIO.to_dict(SeqIO.parse(args.f.split(',')[0],"fasta"))
    second_seqs = SeqIO.to_dict(SeqIO.parse(args.f.split(',')[1],"fasta"))

    first_map = build_sequence_map(first_seqs.keys())
    second_map = build_sequence_map(second_seqs.keys())

    for key in first_map:
        if key in second_map: 

            cur_key_dir = "{}/{}".format(args.o,key)

            make_directory(cur_key_dir)

            # one of these lists should be of size one
            for entry1_id in first_map[key]:
                for entry2_id in second_map[key]:

                    entry1 = first_seqs[entry1_id] # get the Seq object
                    entry2 = second_seqs[entry2_id]

                    entry1_file = "{}/{}.fsa".format(cur_key_dir,entry1.id)
                    entry2_file = "{}/{}.fsa".format(cur_key_dir,entry2.id)

                    if not os.path.isfile(entry1_file):
                        write_fasta(entry1_file,entry1.id,str(entry1.seq))
                    if not os.path.isfile(entry2_file):
                        write_fasta(entry2_file,entry2.id,str(entry2.seq))

                    run_needle(
                        args.n,
                        entry1_file,
                        entry2_file,
                        "{}/{}_WITH_{}.align.txt".format(cur_key_dir,entry1.id,entry2.id)
                    )

                    os.remove(entry1_file) # sequences already stored in inputs
                    os.remove(entry2_file)


def build_sequence_map(sequence_ids):
    """ 
    Builds a map of a base sequence to all of its assembled sequences. For 
    example, ABC.123.1 and ABC.123.2 will be mapped by 123 to those two 
    distinct IDs. 

    Args:
        sequence_id_dict (dict_keys): .keys() from SeqIO.to_dict result

    Returns:
        A dict which maps a base ID like 123 to ABC.123.1 and ABC.123.2. All 
        values in this dict will be a list, most will likely be of length 1.
    """
    sequence_map = defaultdict(list)

    for id in sequence_ids:
        sequence_map[id.split('.')[1]].append(id)

    return sequence_map


def run_needle(needle_exe,aseq,bseq,outfile):
    """ 
    Executes the EMBOSS needle program.

    Args:
        needle_exe (str): path to EMBOSS needle executable
        aseq (str): path to first sequence file
        bseq (str): path to second sequence file
        outfile (str): path to the output file to generate
    """
    needle_cline = NeedleCommandline(
        needle_exe,
        asequence=aseq,
        bsequence=bseq,
        gapopen=10,
        gapextend=0.5,
        outfile=outfile
    )

    stdout,stderr = needle_cline()

    return None


if __name__ == '__main__':
    main()

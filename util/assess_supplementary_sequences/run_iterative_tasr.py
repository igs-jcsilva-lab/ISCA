#!/usr/bin/env python3

"""
A script which runs TASR iteratively on individual sequences as opposed to a 
whole set at once. 

TASR: https://github.com/warrenlr/TASR 

Input: 
    1. Path to the TASR executable. 
    2. A file which has, on each line, the path to a set of reads to assemble. 
    3. A file which has a set of sequences to target assembly for.
    4. A location to a working directory. 

Output: 
    1. Numerous directories (one for each sequence specified in 3.) which has
    all the contents of a typical output from TASR. 

Usage:
    run_iterative_tasr.py -t /path/to/TASR -r /path/to/reads.txt -s /path/to/sequences.fsa -o /path/to/outdir

Author: 
    James Matsumura
"""

import argparse
import os
import subprocess
from Bio import SeqIO
from shared_fxns import make_directory,write_fasta
from shutil import copyfile


def main():

    parser = argparse.ArgumentParser(description='Script to run TASR in an iterative fashion on individual sequences.')
    parser.add_argument('-t', type=str, required=True, help='Location of the TASR exe (~/tasr_v1.6.2/TASR).')
    parser.add_argument('-r', type=str, required=True, help='Path to a file with paths to reads to be used in assembly on each line.')
    parser.add_argument('-s', type=str, required=True, help='Path to a file with the sequences which should be the targets for assembly.')
    parser.add_argument('-o', type=str, required=True, help='Location to generate output directories.')
    args = parser.parse_args()

    make_directory(args.o)

    seqs = SeqIO.to_dict(SeqIO.parse(args.s,"fasta"))

    for id in seqs:

        cur_dir = "{}/{}".format(args.o,id)
        cur_fsa = "{}/seq.fasta".format(cur_dir)
        cur_rds = "{}/reads.txt".format(cur_dir)

        make_directory(cur_dir)
        write_fasta(cur_fsa,id,str(seqs[id].seq))
        copyfile(args.r, cur_rds)

        command = (
            "{} -s {} -f {} -w 1 -u 1 -c 1"
            .format(args.t, cur_fsa, cur_rds)
        )
        
        subprocess.call(command.split())


if __name__ == '__main__':
    main()
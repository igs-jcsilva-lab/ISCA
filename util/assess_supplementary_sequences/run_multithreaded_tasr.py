#!/usr/bin/env python3

"""
A script which runs TASR in a multithreaded environment with each thread 
responsible for assembling an individual sequence.

TASR: https://github.com/warrenlr/TASR 

Input: 
    1. Path to the TASR executable. 
    2. Number of threads to use (this will invoke 1 less than requested here
    to accommodate this running script itself)
    3. A file which has, on each line, the path to a set of reads to assemble. 
    4. A file which has a set of sequences to target assembly for.
    5. A location to a working directory. 

Output: 
    1. Numerous directories (one for each sequence specified in 3.) which each
    haev all the contents of a typical output from TASR. The file of interest
    is likely *contigs in each directory.

Usage:
    run_multithreaded_tasr.py -t /path/to/TASR -d 10 -r /path/to/reads.txt -s /path/to/sequences.fsa -o /path/to/outdir

Author: 
    James Matsumura
"""

import argparse
import multiprocessing as mp
import subprocess
from Bio import SeqIO
from shared_fxns import make_directory,write_fasta
from shutil import copyfile


def main():

    parser = argparse.ArgumentParser(description='Script to run TASR in a multithreaded fashion on individual sequences.')
    parser.add_argument('-t', type=str, required=True, help='Location of the TASR exe (~/tasr_v1.6.2/TASR).')
    parser.add_argument('-d', type=int, required=True, help='Number of threads to use.')
    parser.add_argument('-r', type=str, required=True, help='Path to a file with paths to reads to be used in assembly on each line.')
    parser.add_argument('-s', type=str, required=True, help='Path to a file with the sequences which should be the targets for assembly.')
    parser.add_argument('-o', type=str, required=True, help='Location to generate output directories.')
    args = parser.parse_args()

    make_directory(args.o)

    seqs = SeqIO.to_dict(SeqIO.parse(args.s,"fasta"))

    manager = mp.Manager()
    pool = mp.Pool(args.d-1)
    jobs = []

    for id in seqs:

        cur_dir = "{}/{}".format(args.o,id)
        cur_fsa = "{}/seq.fasta".format(cur_dir)
        cur_rds = "{}/reads.txt".format(cur_dir)

        make_directory(cur_dir)
        write_fasta(cur_fsa,id,str(seqs[id].seq))
        copyfile(args.r, cur_rds)

        jobs.append(pool.apply_async(run_tasr, (args.t,cur_fsa,cur_rds)))

    for job in jobs:  # Get all the returns from the apply_async function
        job.get()
    
    pool.close() #  Tell the queue it's done getting new jobs
    pool.join() # Make sure these new jobs are all finished


def run_tasr(tasr_exe, seqs, reads):
    """ 
    Runs TASR. 

    Args:
        tasr_exe (str): path to the TASR executable
        seqs (str): path to the target sequence file
        reads (str): path to the file with reads to be assembled
    """    
    command = (
        "{} -s {} -f {} -w 1 -u 1 -c 1"
        .format(tasr_exe, seqs, reads)
    )

    subprocess.call(command.split())


if __name__ == '__main__':
    main()
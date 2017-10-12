#!/usr/bin/env python3

"""
This script aims to reformat all directories created by fastq_reads_to_fastq_alleles.py
into something that SGE grid can use for a job array submission. This script
is meant to be run before the grid job to perform assembly is initiated.
The output is a basic map file with two columns, ref and grid ID. The grid IDs will
be arbitrary and start at 1 climbing until it hits the number of directories made.

    Input:
        1. Path to *_ref_map.tsv output from analyze_bam.py
        2. Path to where the output directory for the FASTQs went, same as what 
        was used for fastq_reads_to_fastq_alleles.py
        3. Path to the the directory to initialize directories for all the 
        assembly output
        4. Name of an output map file (maps the ref to the SGE ID))

    Output:
        1. Renamed output directories (locus ID --> random int for grid job)
        2. Output file to know the contents of the renamed directories

    Usage:
        format_for_assembly.py --ref_map /path/to/ref_map.tsv --reads_dir /path/to/reads_dir --assmb_path /path/to/assmb/dirs --outfile /path/to/out_map.tsv

    Author: 
        James Matsumura
"""

import argparse,os,errno
from shared_fxns import make_directory

def main():

    parser = argparse.ArgumentParser(description='Script to set up for SPAdes alignment on a grid.')
    parser.add_argument('--ref_map', '-rm', type=str, required=True, help='Path to *_ref_map.tsv output from analyze_bam.py.')
    parser.add_argument('--reads_dir', '-rd', type=str, required=True, help='Path to where the output directory for the FASTQs went, same as what was used for fastq_reads_to_fastq_alleles.py.')
    parser.add_argument('--assmb_path', '-ap', type=str, required=True, help='Path to the the directory to initialize directories for all the assembly output.')
    parser.add_argument('--outfile', '-o', type=str, required=True, help='Path to output map (maps the ref to the SGE ID)).')
    args = parser.parse_args()

    ref_map = {} # dict to hold the ref and its arbitrary ID starting at 1
    id = 1

    # First, rename the directories
    with open(args.ref_map,'r') as infile:
        for line in infile:
            
            line = line.rstrip()
            ref = line.split('\t') # really just want the first column which is the ref ID

            # Now rename the original directory so that it can be iterated over in a grid
            # job. 
            old_dir = "{0}/{1}".format(args.reads_dir,ref[0])
            new_dir = "{0}/{1}".format(args.reads_dir,id)

            try:
                os.rename(old_dir,new_dir)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            else:
                ref_map[ref[0]] = id

                # Make a new output directory for all the SPAdes assembly files
                spades_out_dir = "{0}/{1}".format(args.assmb_path,id)
                make_directory(spades_out_dir)

                id += 1 # if no exception, it was renamed and need a new ID

    # Now generate a map to know which directories correlate to what IDs
    with open(args.outfile,'w') as outfile:
        for k,v in ref_map.items():
            outfile.write("{0}\t{1}\n".format(k,v))


if __name__ == '__main__':
    main()
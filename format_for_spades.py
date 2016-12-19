

# This script aims to reformat all directories created by fastq_reads_to_fastq_alleles.py
# into something that SGE grid can use for a job array submission. This script
# is meant to be run both before and after the grid job to run SPAdes is compelted.
# The output is a basic map file with two columns, ref and grid ID. The grid IDs will
# be arbitrary and start at 1 climbing until it hits the number of directories made.
#
# Run the script using a command like this:
# python3 analyze_bam.py -i /path/to/ref_map.tsv -path /path/to/ref_dirs -spades_path /path/to/spades/dirs -out /path/to/out_map.tsv
#
# Author: James Matsumura

import argparse,os,errno
from shared_fxns import make_directory

def main():

    parser = argparse.ArgumentParser(description='Script to set up for SPAdes alignment on a grid.')
    parser.add_argument('-ref_map', type=str, required=True, help='Path to *_ref_map.tsv output from analyze_bam.py.')
    parser.add_argument('-path', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-spades_path', type=str, required=True, help='Path to the the directory to initialize directories for all the SPAdes output.')
    parser.add_argument('-outfile', type=str, required=True, help='Path to output map (maps the ref to the SGE ID)).')
    args = parser.parse_args()

    ref_map = {} # dict to hold the ref and its arbitrary ID starting at 1
    id = 1

    # First, rename the directories
    with open(args.refs,'r') as infile:
        for line in infile:
            
            line = line.rstrip()
            ref = line.split('\t') # really just want the first column which is the ref ID

            # Now rename the original directory so that it can be iterated over in a grid
            # job. 
            old_dir = "{0}/{1}".format(args.path,ref[0])
            new_dir = "{0}/{1}".format(args.path,id)

            try:
                os.rename(old_dir,new_dir)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            else:
                ref_map[ref[0]] = id

                # Make a new output directory for all the SPAdes assembly files
                spades_out_dir = "{0}/{1}".format(args.spades_path,id)
                make_directory(spades_out_dir)

                id += 1 # if no exception, it was renamed and need a new ID

    # Now generate a map to know which directories correlate to what IDs
    with open(args.outfile,'w') as outfile:
        for k,v in ref_map.items():
            outfile.write("{0}\t{1}\n".format(k,v))


if __name__ == '__main__':
    main()
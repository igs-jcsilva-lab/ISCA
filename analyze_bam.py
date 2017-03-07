

# This script will parse through a BAM file resulting from a paired-end 
# alignment and extract all query sequence IDs that meet a certain percent
# identity threshold (which the user can set). 
#
# The output will be two tab-delimited files:
# 1) prefix_read_map.tsv
# A map file where a given read is a key and it maps to all the references
# that it hit. The first column is the read ID, every subsequent column is a
# combination of %ID|length|referenceID for all the reference alignments 
# that passed the threshold. Note that %ID and length are included so that
# these values can help determine which read may have a better alignment
# given a similar %ID. 
#
# 2) prefix_ref_map.tsv
# Similar to the first, except this one uses the reference ID as the main
# map point. The subsequent columns will consist of %ID|readID. 
# 
# While these two files are similar, they're more useful in separate cases. The
# first is what will be used downstream here while the second can be used to 
# check for which reference contigs were able to recruit the most reads. 
#
# Run the script using a command like this:
# python3 analyze_bam.py -bam /path/to/bowtie_out.bam (-threshold 90) -o /path/to/out_prefix
#
# Author: James Matsumura

import argparse,pysam
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to isolate all reads and where they aligned to given a BAM file.')
    parser.add_argument('-bam', type=str, help='Path to a BAM file derived from Bowtie2/GSNAP.')
    parser.add_argument('-sam', type=str, help='Path to a SAM file derived from Bowtie2/GSNAP.')
    parser.add_argument('-threshold', type=str, required=False, help='Minimum %ID threshold to retain (entering 95 means %95 minimum %ID). Defaults to %80.')
    parser.add_argument('-o', type=str, required=True, help='Path to prefix of where the two output TSV files should go.')
    args = parser.parse_args()

    if args.sam is None and args.bam is None:
        parser.error("Either -sam or -bam input is required.")
    
    i = None # will be the input file

    if args.sam is not None:
        i = pysam.AlignmentFile(args.sam,'r')
    if args.bam is not None:
        i = pysam.AlignmentFile(args.bam,'rb')

    rd_map,rf_map = (defaultdict(list) for j in range(2)) # establish these dicts as holding lists

    threshold = 80  # minimum cutoff to include a read
    if args.threshold: # user specificying threshold
        threshold = float(args.threshold)

    for read in i.fetch(until_eof=True): # iterate over all reads in the [S|B]AM file.
        que,ref = ("" for j in range(2))

        cigar = read.cigartuples # extract a tuple for the CIGAR string

        # GSNAP may arbitrarily write an invalid start position to denote an invalid
        # alignment. Seems to be particular to the version. 
        if read.reference_start < 0:
            continue

        if cigar is not None: # only act if a CIGAR string is present
            length = read.query_alignment_length # grab length of alignment
            percent_id = calculate_percent_id(cigar,length)
            
            if percent_id >= threshold: # good enough %ID? go on...
                rea = read.query_name # grab read ID
                ref = read.reference_name # grab ID of what the ID aligned to
                rd_map[rea].append("{0}|{1}|{2}".format(percent_id,length,ref))

                # For the reference map, consolidate all GMAP values to a single locus,
                # e.g. make NF54.1.1 and 3D7.1.1 just 1.1
                ele = ref.split('.')
                ref = ele[1]
                rf_map[ref].append("{0}|{1}|{2}".format(percent_id,length,rea))

    i.close()

    # Time to generate the output
    read_map = args.o + "_read_map.tsv"
    ref_map = args.o + "_ref_map.tsv"

    with open(read_map,'w') as o1:
        write_tsv(rd_map,o1)

    with open(ref_map,'w') as o2:
        write_tsv(rf_map,o2)


# Calculates % identity given info from a CIGAR string
# Expects a tuple based on the pysam documentation (http://pysam.readthedocs.io/en/latest/api.html) where...
# 0 is M
# 1 is I
# 2 is D
# 3 is N 
# 4 is S 
# 5 is H 
# 6 is P 
# 7 is =
# 8 is X
# Thus, need to aggregate all 0/M and then divide by total alignment length.
# Arguments:
# cigar = CIGAR tuple
# length = length of the alignment (needed to calculate %ID)
def calculate_percent_id(cigar,length):
    
    matches,percent_id = (0 for j in range(2))
    
    for x in cigar: # iterate over each tuple found
        if x[0] == 0:
            matches += x[1] # only increment match value when 0/M found
    
    percent_id = (100 * (matches/length))

    return float("{0:.2f}".format(percent_id))

# Function which takes the built dictionaries and writes the TSVs
# Arguments:
# map = dictionary of lists where the key is either the read/ref and the list
# contains all the ref/reads that have alignment to the key
# outfile = path to the TSV outfile
def write_tsv(map,outfile):
    for entry in map:
        mapped_vals = "\t".join(map[entry])
        tsv_line = "{0}\t{1}\n".format(entry,mapped_vals)
        outfile.write(tsv_line)

if __name__ == '__main__':
    main()
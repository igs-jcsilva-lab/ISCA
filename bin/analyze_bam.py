

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
# python3 analyze_bam.py --bam /path/to/bowtie_out.bam (-threshold 90) -p /path/to/out_prefix --ea_map /path/to/out_from_extract_alleles.tsv 
#
# Author: James Matsumura

import argparse,pysam
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to isolate all reads and where they aligned to given a BAM file.')
    parser.add_argument('--bam', '-b', type=str, help='Path to a BAM file.')
    parser.add_argument('--sam', '-s', type=str, help='Path to a SAM file.')
    parser.add_argument('--ea_map', '-eam', type=str, help='Path to map.tsv output from extract_alleles.py.')
    parser.add_argument('--threshold', '-t', type=int, default=80, required=False, help='Minimum %ID threshold to retain (entering 95 means %95 minimum %ID). Defaults to %80.')
    parser.add_argument('--prefix', '-p', type=str, required=True, help='Name of the prefix of where the two output TSV files should go.')
    args = parser.parse_args()

    if args.sam is None and args.bam is None:
        parser.error("Either --sam or --bam input is required.")
    
    i = None # will be the input file

    if args.sam is not None: # separate file handlers for noncompressed/compressed
        i = pysam.AlignmentFile(args.sam,'r')
    if args.bam is not None:
        i = pysam.AlignmentFile(args.bam,'rb')

    # Collect the lengths from each reference gene
    ref_lengths = defaultdict(list)
    with open(args.ea_map,'r') as loc_map:

        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')

            for j in range(1,len(ele)): # only care about alleles info
                allele_info = ele[j].split('|')
                length = (int(allele_info[2])-int(allele_info[1])+1) # include starting base
                ref_lengths[allele_info[4]] = length

    rd_map,rf_map = (defaultdict(list) for j in range(2)) # establish these dicts as holding lists

    for read in i.fetch(until_eof=True): # iterate over all reads in the [S|B]AM file.

        cigar = read.cigartuples # extract a tuple for the CIGAR string

        # GSNAP may arbitrarily write an invalid start position to denote an invalid
        # alignment. Seems to be particular to the version. 
        if read.reference_start < 0:
            continue

        if cigar is not None: # only act if a CIGAR string is present

            rea = read.query_name # grab read ID
            ref = read.reference_name # grab ID of what the ID aligned to

            percent_id = calculate_percent_id(cigar,ref_lengths[ref])
  
            if percent_id >= float(args.threshold): # good enough %ID? go on...
                
                # Write out the number of match alignment
                rd_map[rea].append("{0}|{1}|{2}".format(percent_id,read.query_alignment_length,ref))

                # For the reference map, consolidate all GMAP values to a single locus,
                # e.g. make NF54.1.1 and 3D7.1.1 just 1.1
                ele = ref.split('.')
                ref = ele[1]
                rf_map[ref].append("{0}|{1}|{2}".format(percent_id,read.query_alignment_length,rea))

    i.close()

    # Time to generate the output
    read_map = args.prefix + "_read_map.tsv"
    ref_map = args.prefix + "_ref_map.tsv"

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
    
    matches,percent_id,total = (0 for j in range(3))
    
    # If we are using the read length, establish this length dynamically using 
    # all the counts found in the CIGAR
    for x in cigar: # iterate over each tuple found
        if x[0] == 0:
            matches += x[1] # only increment match value when 0/M found
            total += x[1]
        else:
            total += x[1]

    if length < total:
        percent_id = (100 * (matches/length))
    else:
        percent_id = (100 * (matches/total))

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
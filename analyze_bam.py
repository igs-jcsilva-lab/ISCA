

# This script will parse through a BAM file resulting from a paired-end 
# alignment and extract all query sequence IDs that meet a certain percent
# identity threshold (which the user can set). The output will be a tab-delimited 
# file with the first column being the read ID, the second column being the
# %ID found (calculated by CIGAR string in SAM/BAM file), and all subsequent
# columns being the reference that the read aligned to.  
#
# Run the script using a command like this:
# python3 analyze_bam.py -bam /path/to/bowtie_out.bam (-threshold 90) -o /path/to/output.tsv
#
# Author: James Matsumura

import sys,re,argparse,pysam

def main():

    parser = argparse.ArgumentParser(description='Script to isolate all reads and where they aligned to given a BAM file.')
    parser.add_argument('-bam', type=str, required=True, help='Path to a BAM file derived from Bowtie2.')
    parser.add_argument('-threshold', type=str, required=False, help='Minimum %ID threshold to retain (entering 95 means %95 minimum %ID). Defaults to %80.')
    parser.add_argument('-o', type=str, required=True, help='Path to where the output TSV should go.')
    args = parser.parse_args()

    i = pysam.AlignmentFile(args.bam,'r')
    o = open(args.o,'w')

    threshold = 80  # minimum cutoff to include a read
    if args.threshold: # user specificying threshold
        threshold = float(args.threshold)

    for read in i: # iterate over all reads in the [S|B]AM file.
        que,ref = ("" for i in range(2))

        cigar = read.cigartuples # extract a tuple for the CIGAR string

        if cigar is not None: # only act if a CIGAR string is present
            length = read.query_alignment_length # grab length of alignment
            percent_id = calculate_percent_id(cigar,length)
            
            if percent_id >= threshold: # good enough %ID? go on...
                que = read.query_name # grab read ID
                ref = read.reference_name # grab ID of what the ID aligned to


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
def calculate_percent_id(cigar,length):
    
    matches,percent_id = (0 for i in range(2))
    
    for x in cigar: # iterate over each tuple found
        if x[0] == 0:
            matches += x[1] # only increment match value when 0/M found
    
    percent_id = (100 * (matches/length))

    return float("{0:.2f}".format(percent_id))

if __name__ == '__main__':
    main()
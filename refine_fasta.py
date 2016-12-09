

# This follows analyze_bam.py and accepts the *_read_map.tsv generated from that,
# NOT the *_ref_map.tsv. Remember that this dataset must consist of PAIRED-END 
# reads. This script expects the reads to match in their name except for the final
# number value, e.g. ABC.123.1 + ABC.123.2 are a properly formatted pair.
#
# This script aims to accomplish two things:
#
# 1) Generate some statistics on the alignments like how many reads mapped 
# to multiple loci, how many reads mapped to only aligned to one locus, and
# how many reads have a discrepancy where the two mates in a pair map to 
# different loci. All these stats will be written to STDOUT.
#
# 2) Using the alignment information, potentially refine a FASTA file of reads
# to either those pairs found aligning to just one locus, or no filtering at all.
#
# Run the script using a command like this:
# python3 analyze_bam.py -i /path/to/analyze_bam.out -f /path/to/reads.fsa -filter (yes|no) -o /path/to/out.fsa
#
# Author: James Matsumura

import sys,re,argparse
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to generate stats given output from analyze_bam.py and filter a set of paired-end FASTA reads.')
    parser.add_argument('-i', type=str, required=True, help='Path to *_read_map.tsv output from analyze_bam.py.')
    parser.add_argument('-f', type=str, required=True, help='Path to the original FASTQ file of paired-end reads.')
    parser.add_argument('-filter', type=str, required=True, help='Either "yes" or "no" for removing discrepancies + multi-locus mapping reads.')
    parser.add_argument('-o', type=str, required=True, help='Path to where the output FASTA should go.')
    args = parser.parse_args()

    i = open(args.i,'r')
    f = open(args.f,'r')
    filter = args.filter
    o = open(args.o,'w')

    counts = {'single_map':0,'multi_map':0,'discrepancy':0} # count these stats as they are processed. 

    r1,r2 = (defaultdict(list) for i in range(2)) # establish each mate dict as an empty list
    # Establish two sets...
    # One to keep track of which mates have been checked to make sure none are missed
    # One to keep track of the shared ID across mates for filtering the FASTQ input
    checked_ids,ids_to_keep = (set() for i in range(2)) 

    # This first iteration only cares about grabbing all mates and their reference alignment stats
    for line in i: 
        
        line = line.rstrip()
        ele = line.split('\t')
        
        if ele[0][-1] == "1": # read mate 1
            for x in range(1,len(ele)):
                r1[ele[0]].append(ele[x])
        else: # read mate 2
            for x in range(1,len(ele)):
                r2[ele[0]].append(ele[x])

    shared_id = "" # id in the format of ABC.123 for pairs ABC.123.1 + ABC.123.2

    # Now, iterate over each dict of mates and filter if required
    for read in r1: # mate 1
        shared_id = read[:-2]
        mate_id = shared_id + ".2"
        count_val = verify_alignment(r1[read],r2[mate_id])
        counts[count_val] += 1

        if filter == "yes": # need to isolate best allele across pairs
            if count_val == "single_map":
                ids_to_keep.add(shared_id)

        checked_ids.add(shared_id) # identify these as looked at before going into r2 dict

    for read in r2: # mate 2, only check if the mate wasn't caught by the r1 dict
        shared_id = read[:-2]

        if shared_id not in checked_ids: # if not checked using mate 1, verify now.

            mate_id = shared_id + ".1"   
            count_val = verify_alignment(r1[mate_id],r2[read])
            counts[count_val] += 1

            if filter == "yes":    
                if count_val == "single_map":
                    ids_to_keep.add(shared_id)

            checked_ids.add(shared_id)

    if filter == "no": # if no filtering occurred, establish ids_to_keep as equivalent to all checked
        ids_to_keep = checked_ids

    for k,v in counts: # give the user some idea of how much they are potentially filtering out
        print("{0} read-pairs have a {1}.\n".format(k,v))

# Function to compare where the two mates in a pair mapped to. Returns 
# 'single_map' if both only map to a single locus, 'multi_map' if one
# or both of the reads map to more than one locus, and 'discrep'  if 
# the two mates do not map to the same locus. 
# Arguments:
# list1 = list of alignments from the first mate
# list2 = list of alignments from the second mate
def verify_alignment(list1,list2):

    set1 = isolate_loci(list1)
    set2 = isolate_loci(list2)
    
    if len(set1) == 1 and set1 == set2: # only one locus, and the same one in both mates
        return "single_map"
    elif (len(set1) == 1 and set2 is None): # this and the next account for where one read maps and the other doesn't
        return "single_map"
    elif (len(set2) == 1 and set1 is None):
        return "single_map"
    elif set1 != set2: # not sharing the same loci, discrepancy!
        return "discrepancy"
    elif (len(set1) > 1 or len(set2 > 1)) and set1 == set2: # both reads mapping to more than one locus
        return "multi_map" 

# Function which will return a set of all unique loci found given a list of
# alignment info.
# Arguments:
# alignment_list = expects the following data in each element: %ID|length|reference
def isolate_loci(alignment_list):

    loci = set() # set of reference IDs

    for ref in alignment_list:
        
        ele = ref.split('|')
        align_info = ele[2].split('.') # grab just the reference portion
        id = align_info[1] # second element is the reference without pre/suf-fixes

        loci.add(id)

    return loci


if __name__ == '__main__':
    main()
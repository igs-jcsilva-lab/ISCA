

# This follows analyze_bam.py and accepts the *_read_map.tsv generated from that 
# script. Remember that this dataset must consist of PAIRED-END SORTED FASTQ reads
# that have the SAME order. This is because this script will iterate over the two
# files simultaneously. This script expects the reads to match in their name except for 
# the final number value, e.g. ABC.123.1 + ABC.123.2 are a properly formatted pair.
#
# This script aims to accomplish two things:
#
# 1) Generate some statistics on the alignments like how many reads mapped 
# to multiple loci, how many reads mapped to only aligned to one locus, and
# how many reads have a discrepancy where the two mates in a pair map to 
# different loci. All these stats will be written to STDOUT.
#
# 2) Using the alignment information, potentially refine a FASTQ file of reads
# to either those pairs found aligning to just one locus, or no filtering at all.
#
# The output will be a set of directories, one for each locus aligned to by at least
# one read. Under each directory is the set of paired reads that an alignment
# originated from.
#
# Run the script using a command like this:
# python3 analyze_bam.py -i /path/to/analyze_bam.out -f /path/to/reads.fsa -filter (yes|no) -o /path/to/out.fsa
#
# Author: James Matsumura

import sys,re,argparse,gzip,os,itertools
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to generate stats given output from analyze_bam.py and filter a set of paired-end FASTQ reads.')
    parser.add_argument('-reads', type=str, required=True, help='Path to *_read_map.tsv output from analyze_bam.py.')
    parser.add_argument('-fastq', type=str, required=True, help='Path to the original FASTQ file prefix of paired-end reads (e.g., enter ABC.123 for pairs ABC.123.1+ABC.123.2). MUST be gunzipped.')
    parser.add_argument('-filter', type=str, required=True, help='Either "yes" or "no" for removing discrepancies + multi-locus mapping reads.')
    parser.add_argument('-out', type=str, required=True, help='Path to where the output directory for the FASTQs to go.')
    args = parser.parse_args()

    filename = args.fastq 
    filter = args.filter
    output = args.out

    counts = {'single_map':0,'multi_map':0,'discrepancy':0} # count these stats as they are processed. 

    # Establish three dicts:
    # first two dicts consist of one for each mate
    # third dict is the IDs that need to be mapped (checking based on if the user wants to filter)
    r1,r2,ids_to_keep = (defaultdict(list) for j in range(3)) # establish each mate dict as an empty list

    # This first iteration only cares about grabbing all mates and their reference alignment info
    with open(args.reads,'r') as reads:
        for line in reads: 
            
            line = line.rstrip()
            ele = line.split('\t')
            
            if ele[0][-1] == "1": # read mate 1
                for j in range(1,len(ele)):
                    alignment = ele[j].split('|') # split the alignment data
                    ref = alignment[2].split('.') # split the reference name
                    ref_loc = ref[1] # grab just the base reference locus
                    # don't double up on references (possible if mapping to same locus from different samples)
                    if ref_loc not in r1[ele[0]]: 
                        r1[ele[0]].append(ref_loc)
            else: # read mate 2
                for j in range(1,len(ele)):
                    alignment = ele[j].split('|')
                    ref = alignment[2].split('.')
                    ref_loc = ref[1]
                    if ref_loc not in r2[ele[0]]:
                        r2[ele[0]].append(ref_loc)

    shared_id = "" # id in the format of ABC.123 for pairs ABC.123.1 + ABC.123.2
    checked_ids = set() # set to speed up processing of R2 if already covered by R1

    # Now, iterate over each dict of mates and filter if required
    for read in r1: # mate 1
        shared_id = read[:-2]
        mate_id = shared_id + ".2"

        # Generate stats regardless of filtering or not, can help the user decide if they should
        count_val = verify_alignment(r1[read],r2[mate_id])
        counts[count_val] += 1

        if filter == "yes" and count_val == "single_map": # need to isolate reads that only map once
            
            # If a single map value, know that both reads share the same locus
            if not r1[read]: # if R1 didn't map, means R2 did
                ids_to_keep[shared_id].append(r2[mate_id])
            elif not r2[mate_id]: # same as above, if R2 didn't map, means R1 did
                ids_to_keep[shared_id].append(r1[read])
            else: # else, they both mapped to the same locus and can use either value
                ids_to_keep[shared_id].append(r1[read])

        else: # no filter needed, add all distinct loci found per read

            for ref in r1[read]:
                ids_to_keep[shared_id].append(ref)
            for ref in r2[mate_id]:
                if ref not in ids_to_keep[shared_id]: # make sure not to double up on loci across mates
                    ids_to_keep[shared_id].append(ref)

        checked_ids.add(shared_id) # identify these as looked at before going into r2 dict

    for read in r2: # mate 2, only check if the mate wasn't caught by the r1 dict
        shared_id = read[:-2]

        if shared_id not in checked_ids: # if not checked using mate 1, verify now.

            mate_id = shared_id + ".1"   
            count_val = verify_alignment(r1[mate_id],r2[read])
            counts[count_val] += 1

            # If we are here, the read was not found in R1. Thus, get loci strictly from R2.
            if filter == "yes" and count_val == "single_map": 
                ids_to_keep[shared_id].append(r2[read])

            else: # Again, was not found in R1 so we know all loci are from R2. 
                for ref in r2[read]:
                    ids_to_keep[shared_id].append(ref)

    # At this point, ids_to_keep now has a dictionary mapping all read IDs to loci that
    # they aligned to. This is all that's needed to build a set of directories that house
    # reads just mapping to those loci for use in assembly.

    r1,r2,checked_ids = (None for j in range(3)) # done with these, free up some memory

    # give the user some idea of how much they are potentially filtering out
    out_stats = output + ".stats"
    with open(out_stats,'w') as stats_file:
        for k,v in counts.items():
            stats_file.write("{0} read-pairs have a {1}.\n".format(v,k))

    # Regardless of filtering based on alignment single/multiple/discrepancies or not, still
    # need to filter all the FASTQ reads to just those that aligned to a gene region.
    filter_fastq(ids_to_keep,filename,output)


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
    
    # only one locus, and the same one in both mates
    if len(set1) == 1 and set1 == set2: 
        return "single_map"
    # this and the next account for where one read maps and the other doesn't
    elif (len(set1) == 1 and set2 is None): 
        return "single_map"
    elif (len(set2) == 1 and set1 is None):
        return "single_map"
    # some reads are mapping to more than one locus
    elif ((len(set1) > 1) or (len(set2) > 1)): 
        return "multi_map" 
    # simmply not sharing the same loci, discrepancy!
    elif set1 != set2: 
        return "discrepancy"

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

# Function to parse through a FASTQ file and generate new ones that only consist
# of IDs, per locus, found to be valid by the alignment and this script.
# Arguments:
# ids = set of IDs to be checked against while parsing the FASTQ file.
# fastq = prefix for the input FASTQ files. 
# outdir = directory prefix for where the output will be written. 
def filter_fastq(ids,fastq,outdir):

    file1 = fastq + "1.fastq.gz"
    file2 = fastq + "2.fastq.gz"
    entry1,entry2 = ([] for j in range(2))
    lineno = 0
    seen = 0 # count the reads to potentially leave the files early if all found
    total = len(ids)

    # Iterate over each file simultaneously
    with gzip.open(file1,'rt') as f1:
        with gzip.open(file2,'rt') as f2:
            for line1,line2 in zip(f1,f2):
                entry1.append(line1)
                entry2.append(line2)
                lineno += 1

                if lineno == 4: # got a FASTQ entry, check if it's relevant
                    # Note that these mates will both be included if just one is relevant,
                    # so can do all checks using just one of the mates.
                    header = entry1[0]
                    elements = header.split(' ')
                    id = elements[0][1:] # drop the '@'
                    id = id[:-2] # drop the mate distinction of '.1' or '.2'

                    if id in ids: # if relevant, write to all necessary directories/files
                        seen += 1

                        # Establish all loci mapped to, could be many if not filtering
                        for ref in ids[id]:
                            dir = "{0}/{1}".format(outdir,ref)
                            if not os.path.exists(dir): # create dir if not present already
                                os.makedirs(dir)

                            out1 = dir + "/R1.fastq.gz"
                            out2 = dir + "/R2.fastq.gz"

                            # add to whatever FASTQ file is already there
                            with gzip.open(out1,'ab') as o1:
                                for l in entry1:
                                    o1.write(l.encode())
                            with gzip.open(out2,'ab') as o2:
                                for l in entry2:
                                    o2.write(l.encode())


                    entry1,entry2 = ([] for j in range(2)) # reset for next entry
                    lineno = 0

                if seen == total: # got them all, leave
                    break


if __name__ == '__main__':
    main()
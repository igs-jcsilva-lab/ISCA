

# This follows analyze_bam.py and accepts the *_read_map.tsv and the *_ref_map.tsv
# generated from that script. Remember that this dataset must consist of PAIRED-END 
# FASTQ reads. This script expects the reads to match in their name except for 
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
# The output will be a set of directories, one for each locus identified in the 
# *_ref_map.tsv file and containing the set of paired reads that aligned to it 
# in one of the reference allele sequences. 
#
# Run the script using a command like this:
# python3 analyze_bam.py -i /path/to/analyze_bam.out -f /path/to/reads.fsa -filter (yes|no) -o /path/to/out.fsa
#
# Author: James Matsumura

import sys,re,argparse,gzip,os
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to generate stats given output from analyze_bam.py and filter a set of paired-end FASTQ reads.')
    parser.add_argument('-reads', type=str, required=True, help='Path to *_read_map.tsv output from analyze_bam.py.')
    parser.add_argument('-refs', type=str, required=True, help='Path to *_ref_map.tsv output from analyze_bam.py.')
    parser.add_argument('-fastq', type=str, required=True, help='Path to the original FASTQ file prefix of paired-end reads (e.g., enter ABC.123 for pairs ABC.123.1+ABC.123.2). MUST be gunzipped.')
    parser.add_argument('-filter', type=str, required=True, help='Either "yes" or "no" for removing discrepancies + multi-locus mapping reads.')
    parser.add_argument('-out', type=str, required=True, help='Path to where the output directory for the FASTQs to go.')
    args = parser.parse_args()

    filter = args.filter
    output = args.out

    counts = {'single_map':0,'multi_map':0,'discrepancy':0} # count these stats as they are processed. 

    r1,r2 = (defaultdict(list) for i in range(2)) # establish each mate dict as an empty list
    # Establish two sets...
    # One to keep track of which mates have been checked to make sure none are missed
    # One to keep track of the shared ID across mates for filtering the FASTQ input
    checked_ids,ids_to_keep = (set() for i in range(2)) 

    # This first iteration only cares about grabbing all mates and their reference alignment info
    reads = open(args.reads,'r')

    for line in reads: 
        
        line = line.rstrip()
        ele = line.split('\t')
        
        if ele[0][-1] == "1": # read mate 1
            for x in range(1,len(ele)):
                r1[ele[0]].append(ele[x])
        else: # read mate 2
            for x in range(1,len(ele)):
                r2[ele[0]].append(ele[x])

    reads.close()

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
    
    r1,r2 = (None for i in range(2)) # try free up a bit of memory

    # give the user some idea of how much they are potentially filtering out
    out_stats = output + ".stats"
    stats_file = open(out_stats,'w')
    for k,v in counts.items():
        stats_file.write("{0} read-pairs have a {1}.\n".format(v,k))
    stats_file.close

    # if no filtering occurred, establish ids_to_keep as equivalent to all checked
    if filter == "no": 
        ids_to_keep = checked_ids

    # Regardless of filtering based on alignment single/multiple/discrepancies or not, still
    # need to filter all the FASTQ reads to just those that aligned to a gene region.
    filename = args.fastq 

    # Iterate over the reference map, for each locus isolate the set that is shared between
    # those aligned and within the ids_to_keep set.
    refs = open(args.refs,'r')

    for ref in refs:

        loc_reads = set() # Initialize a new set of reads to keep for each reference

        ele = ref.split('\t')

        # Make a new directory for the locus
        new_directory = "{0}/{1}".format(output,ele[0])
        os.makedirs(new_directory)

        # Iterate over all the reads, for this locus, and check if they're in ids_to_keep
        for i in range(1,len(ele)):
            alignment = ele[i].split('|')
            
            if alignment[2] in ids_to_keep: # if present, get this read for this locus assembly
                loc_reads.add(alignment[2])

        # First, do mate 1
        file1 = filename + "1.fastq.gz" 
        out1 = new_directory + "/R1.fastq.gz"
        f1 = gzip.open(file1,'rt')
        o1 = gzip.open(out1,'wb')

        filter_fastq(loc_reads,f1,o1) 

        f1.close
        o1.close

        # Now mate 2
        file2 = filename + "2.fastq.gz" 
        out2 = new_directory + "/R2.fastq.gz"
        f2 = gzip.open(file2,'rt')
        o2 = gzip.open(out2,'wb')

        filter_fastq(loc_reads,f2,o2) 

        f2.close
        o2.close

    refs.close()


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

# Function to parse through a FASTQ file and generate a new one that only consists
# of IDs found to be valid by the alignment and this script.
# Arguments:
# ids = set of IDs to be checked against while parsing the FASTQ file.
# fastq = handle for the input FASTQ file. 
# outfile = handle for the output FASTQ file. 
def filter_fastq(ids,fastq,outfile):

    lineno = 0 # iterate through each FASTQ entry, 4 lines at a time
    entry = []
    seen = 0 # keep count to potentially leave FASTQ early if seen all IDs
    num_of_reads = len(ids)
    
    for line in fastq:

        entry.append(line)
        lineno += 1

        if lineno == 4: # got an entry, check if it's a relevant one
            header = entry[0]
            elements = header.split(' ')
            id = elements[0][1:] # drop the @
            id = id[:-2] # drop the last two characters which designate which mate

            # Finally, check whether or not ID is in set. Output if so.
            if id in ids:
                seen += 1
                for l in entry:
                    outfile.write(l.encode())

            entry = []
            lineno = 0

        if seen == num_of_reads: # leave early if got all the reads needed
            break


if __name__ == '__main__':
    main()
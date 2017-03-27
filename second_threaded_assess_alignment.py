

# This script parses through the output directories of threaded_global_alignment.py
# to extract the best alignment for each assembled sequence. Note the 'priority' 
# parameter which prefers that isolate set which was used to map from via GMAP. 
# When the desired reference priority is not the best hit, brief stats on what %ID 
# it was capable of reaching for that locus compared to the best hit are presented.
# 
# The STDOUT consists of the following columns (tab-separated):
# locus best_ref %id_of_best_ref priority_ref %id_of_priority_ref %id_difference
#
# The ids_v_cov.tsv file consists of the best hits for each assembled locus. 
# The columns are as follows (tab-separated):
# %ID coverage(reference/assembled) reference path_to_best_alignment
#
# This script differs from first_threaded_assess_alignment.py in that it 
# does some additional analysis on any embedded unaligned regions within
# the assembled sequence compared to the reference sequence. In addition 
# to the columns mentioned above for ids_v_cov.tsv, there will be subsequent 
# columns (if applicable) for any embedded regions which will tell:
# where the regions are and how long the regions are. If there are multiple
# embedded regions then the output will be separated by commas. For instance:
# region1,region2 length1,length2 
#
# Run the script using a command like this:
# python3 second_threaded_assess_alignment.py -assmb_map /path/to/format_for_assembly.tsv -ga_stdout threaded_global_alignment_stdout.txt -algn_path /path/to/alignments -out /path/to/output_dir -priority 3D7
#
# Author: James Matsumura

import re,argparse,os,collections
import multiprocessing as mp
from Bio import AlignIO

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    parser.add_argument('-assmb_map', type=str, required=True, help='Path to map.tsv output from format_for_assembly.py or final_verdict.py.')
    parser.add_argument('-cpus', type=int, required=True, help='Number of cores to use.')
    parser.add_argument('-algn_path', type=str, required=True, help='Path to the the directory preceding all the alignment directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-out', type=str, required=True, help='Path to output directory for these stats.')
    parser.add_argument('-priority', type=str, required=False, help='Optional prefix for prioritizing one isolate over the others.')
    args = parser.parse_args()

    # Set up the multiprocessing manager, pool, and queue
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(args.cpus)
    pool.apply_async(listener, (q,args.out))
    jobs = []

    # If priority is provided, set that value here
    priority = ""
    if args.priority:
        priority = args.priority

    # Need to iterate over the map generated from HGA+SB step.
    with open(args.assmb_map,'r') as loc_map:
        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[0]

            algn_dir = "{0}/{1}".format(args.algn_path,locus)

            job = pool.apply_async(worker, (algn_dir,locus,priority,q))
            jobs.append(job)

    # Get all the returns from the apply_async function.
    for job in jobs:
        job.get()
    
    q.put('stop') # should be no more messages
    pool.close() #  Tell the queue it's done getting new jobs
    pool.join() # Make sure these new jobs are all finished

# This is the worker that each CPU will process asynchronously
# Arguments:
# algn_dir = the locus that HGA+SB attempted to assemble
# locus = particular locus being assessed right now
# priority = if provided, same as args.priority
# queue = queue used to send writes to the outfile
def worker(algn_dir,locus,priority,queue):
    isos,scores,ids,files,cov,nogap_id = ([] for i in range(6)) # reinitialize for every locus

    # If the minimum threshold is set high enough, it is possible for
    # no alignments to have been performed. Print to STDOUT in case
    # this does happen. 
    aligned = False
    
    if not os.path.isdir(algn_dir):
        print("The locus {0} could build a scaffold but failed to find an alignment.".format(locus))
        return

    # Found the alignment directory for this locus, now iterate over 
    # the final alignments and pull the best score.
    for file in os.listdir(algn_dir):
        a,b = ("" for i in range(2)) # store lengths of the trimmed alignments

        if 'Scaffold' in file and file.endswith(".trimmed_align.txt"):

            # If we know which reference we want to assemble, skip all other files. 
            if priority != "" and not file.startswith(priority):
                continue

            aligned = True 
            
            isolate = file.split('.')[0] # grab the reference group
            full_path = "{0}/{1}".format(algn_dir,file)

            # Extract the sequence lengths to establish a ratio of
            # potential coverage. >1 means reference is longer than
            # assembled seq and <1 means the assembled seq is longer.
            alignment = AlignIO.read(full_path,'emboss')
            for sequence in alignment:
                if a == "":
                    a = str(sequence.seq)
                else:
                    b = str(sequence.seq)

                if a != "" and b != "":
                    # Check how many bases of A are covered by B with exact 
                    # matches and output this percentage. Ignore gaps.
                    nogap_id.append(calculate_exact_alignment(a,b)) 

                    # Just get the length of the sequences to calculate coverage.
                    # Note that the presence of spacers or extraneous repeats 
                    # can have a significant impact on shifting the coverage 
                    # ratio to find the assembly as much longer. 
                    a = a.replace('-','')
                    b = b.replace('-','')

            stats = parse_alignment(full_path)

            # Seems getting max from a list is faster than dict
            isos.append(isolate)
            scores.append(float(stats['score']))
            ids.append(stats['id'])
            cov.append(len(a)/len(b))
            files.append(full_path)

    # If no trimmed_align.txt files found, no alignments were performed
    # even though contigs were present.
    if aligned == False:
        print("The locus {0} could build a scaffold but failed to find an alignment.".format(locus))
        return

    best = ids.index(max(ids))
    
    # This block is not needed for the current set of test cases but likely
    # will be needed in the future. 
    # Make sure a tie goes to the prioritized isolate.
    #if priority != "":
    #    m = max(ids)
    #    x = [i for i, j in enumerate(ids) if j == m]
    #    if len(x) > 1: # only a concern in the case of a tie
    #        for i in x:
    #            if isos[i] == priority:
    #                best = i

    best_iso = isos[best]
    best_id = ids[best]
    best_cov = cov[best]
    best_file = files[best]
    best_nogap_id = nogap_id[best]

    queue.put("{0}\t{1}\t{2}\t{3}\t{4}\n".format(best_id,best_cov,best_iso,best_file,best_nogap_id))

    # This block is not needed for the current set of test cases but likely
    # will be needed in the future. 
    #if priority != "":
        # If the best hit was not the prioritized isolate, then print to STDOUT
        # what the best value was for the priority as well as the actual best.
    #    if best_iso != priority:
    #        prioritized_indexes = [i for i,j in enumerate(isos) if j == priority]
    #        prioritized_ids = [ids[i] for i in prioritized_indexes]
            # Make sure the prioritized isolate has an alignment. 
    #        prioritized_best_id = 0
    #        if len(prioritized_ids) > 0:
    #            prioritized_best_id = max(prioritized_ids)
    #        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2f}".format(locus,best_iso,best_id,priority,prioritized_best_id,best_id-prioritized_best_id))

# This will act as the sole writer to the output file. This way there is no 
# concern with locks and what not. 
# Arguments:
# queue = queue used to communicate what should be written out
# out_dir = location of where to write out the file
def listener(queue,out_dir):

    # Listens for messages and writes to the final map file
    outfile = "{0}/second_ids_v_cov.tsv".format(out_dir)
    while 1:
        msg = queue.get()
        if msg == 'stop':
            break
        with open(outfile,'a') as out:
            out.write(str(msg))
            out.flush()

# Function to parse over the output of EMBOSS's Needle program and extract the
# score of the alignment.
# Argument:
# infile = *.trimmed_align.txt file generated from a Needle alignment. 
def parse_alignment(infile):

    stats = {'score':0,'id':0}

    with open(infile,'r') as alignment:
        for line in alignment:
            if line.startswith('# Score:'):
                stats['score'] = re.search(r'#\sScore:\s(.*)$',line).group(1) 
            elif line.startswith('# Identity:'):
                stats['id'] = re.search(r'#\sIdentity:\s+\d+/\d+\s\(\s?(\d+\.\d+)%\)$',line).group(1)
                stats['id'] = float(stats['id'])
            elif line.startswith('a.trimmed'): # reached actual alignment, no need
                break

    return stats

# Function to check how many bases from the reference are mapping to the 
# assembled sequence. 
def calculate_exact_alignment(aseq,bseq):

    total,perfect_match = (0 for i in range(2))

    for a_base,b_base in zip(aseq,bseq):
        # If it's not a gap for A, check if it aligns perfectly to B
        if a_base != "-": # only care about what exists in the reference
            if a_base == b_base: # if match, note it as such
                perfect_match += 1
                total += 1
            else: # not a reference gap and not a match
                total += 1

    return "{0:.2f}".format(perfect_match/total*100)

if __name__ == '__main__':
    main()

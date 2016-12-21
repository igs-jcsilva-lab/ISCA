

# This script follows global_alignment.py and shares some of the same inputs.
# The output will be a stats file describing things like the proportion of
# assembled reads that best mapped to which reference isolates. 
#
# Run the script using a command like this:
# python3 assess_alignment.py -ffs_map /path/to/format_for_spades.tsv -assmb_path -/path/to/spades_out -out /path/to/stats.txt
#
# Author: James Matsumura

import re,argparse,os

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    parser.add_argument('-ffs_map', type=str, required=True, help='Path to map.tsv output from format_for_spades.py.')
    parser.add_argument('-ga_map', type=str, required=False, help='Optional path to map.tsv output from global_alignment.py. If provided, this will filter only by those best alignments where the reference length is entirely captured by the assembly.')
    parser.add_argument('-algn_path', type=str, required=True, help='Path to the the directory preceding all the alignment directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-out', type=str, required=True, help='Path to output directory for these stats.')
    args = parser.parse_args()

    # dict for counting which isolates have the best alignments
    isolate_counts = {} 
    # dict for counting how strong the %ID is for the best alignment
    percent_id = {"0<=x<10":0,"10<=x<20":0,"20<=x<30":0,"30<=x<40":0,
                "40<=x<50":0,"50<=x<60":0,"60<=x<70":0,"70<=x<80":0,
                "80<=x<90":0,"90<=x<100":0,"x=100":0}  
    
    percent_id_list = ["x=100","90<=x<100","80<=x<90","70<=x<80","60<=x<70",
                    "50<=x<60","40<=x<50","30<=x<40","20<=x<30","0<=x<10",
                    "10<=x<20"]

    # A set that, if ga_map is provided, will filter the best alignments where
    # the length of the locus is entirely captured by the asesmbly. 
    locus_assembled = set()

    # First build a set for all the alignments where the locus was entirely
    # captured by the assembly. 
    if args.ga_map:
        with open(args.ga_map,'r') as align_map:
            for line in align_map:
                line = line.rstrip()
                ele = line.split('\t')
                if ele[1] == 'assmb':
                    locus_assembled.add(ele[0])

    # Need to iterate over the map generated from SPAdes step, these refs are
    # guaranteed to have been aligned. 
    with open(args.ffs_map,'r') as loc_map:
        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[0]
            algn_dir = "{0}/{1}".format(args.algn_path,locus)
            isos,scores,ids,filenames = ([] for i in range(3)) # reinitialize for every locus

            # Found the alignment directory for this locus, now iterate over 
            # the final alignments and pull the best score.
            for file in os.listdir(algn_dir):
                if file.endswith(".trimmed_align.txt"):

                    full_path = "{0}/{1}".format(algn_dir,file)

                    isolate = file.split('.')[0] # grab the isolate name
                    stats = parse_alignment(full_path)

                    # Seems getting max from a list is faster than dict
                    isos.append(isolate) 
                    scores.append(float(stats['score']))
                    ids.append(stats['id'])
                    filenames.append(file)

            best = scores.index(max(scores))
            best_iso = isos[best]
            best_id = ids[best]
            best_file = filenames[best]

            # If filtering is necessary, check if this file contains an alignment
            # where the reference was entirely re-assembled by the assembler. 
            if args.ga_map:
                if best_file not in locus_assembled:
                    continue # go to the next set of alignments if not the case

            if best_iso in isolate_counts:
                isolate_counts[best_iso] += 1
            else:
                isolate_counts[best_iso] = 1

            percent_id = bin_percent_id(percent_id,best_id)

    outfile = "{0}/alignment_stats.txt".format(args.out)
    with open(outfile,'w') as out:

        out.write("Number of times a given isolate was best aligned to (if an isolate doesn't appear, never was the best alignment):\n")
        tot = 0
        for k,v in isolate_counts.items():
            tot += v
        for k,v in isolate_counts.items():
            rel = float("{0:.2f}".format(100*v/tot)) # calculate relative percentage
            out.write("{0}\t\t{1} ({2}%)\n".format(k,v,rel))

        out.write("\nNumber of times a given percent identity was found for the best alignment:\n")
        tot = 0
        for k,v in percent_id.items():
            tot += v
        for bin in percent_id_list:
            rel = float("{0:.2f}".format(100*percent_id[bin]/tot))
            if bin != 'x=100' and bin != '0<=x<10':
                out.write("{0}\t\t{1} ({2}%)\n".format(bin,percent_id[bin],rel))
            else:
                out.write("{0}\t\t\t{1} ({2}%)\n".format(bin,percent_id[bin],rel))


# Function to parse over the output of EMBOSS's Needle program and extract the
# score of the alignment.
# Argument:
# infile = *.trimmed_align.txt file generated from a Needle alignment. 
def parse_alignment(infile):

    regex_for_score = r'#\sScore:\s(.*)$'
    regex_for_id = r'#\sIdentity:\s+\d+/\d+\s\(\s?(\d+)\.\d+%\)$'
    stats = {'score':0,'id':0}

    with open(infile,'r') as alignment:
        for line in alignment:
            if line.startswith('# Score:'):
                stats['score'] = re.search(regex_for_score,line).group(1) 
            elif line.startswith('# Identity:'):
                stats['id'] = re.search(regex_for_id,line).group(1)
            elif line.startswith('a.trimmed'): # reached actual alignment, no need
                break

    return stats

# Function to house what is essentially a switch statement for grouping together
# different percent identity values for output.
# Argument:
# id_dict = a percent_id dict that will be built upon
# percent = percent value to use to bin
def bin_percent_id(id_dict,percent):

    percent = int(percent) # can trim decimals due to bin

    # Take advantage of the order of processing to slim down the if statements
    if percent == 100:
        id_dict['x=100'] += 1
    elif percent >= 90:
        id_dict['90<=x<100'] += 1
    elif percent >= 80:
        id_dict['80<=x<90'] += 1
    elif percent >= 70:
        id_dict['70<=x<80"'] += 1
    elif percent >= 60:
        id_dict['60<=x<70'] += 1
    elif percent >= 50:
        id_dict['50<=x<60'] += 1
    elif percent >= 40:
        id_dict['40<=x<50'] += 1
    elif percent >= 30:
        id_dict['30<=x<40'] += 1
    elif percent >= 20:
        id_dict['20<=x<30'] += 1
    elif percent >= 10:
        id_dict['10<=x<20'] += 1
    elif percent >= 0:
        id_dict['0<=x<10'] += 1

    return id_dict


if __name__ == '__main__':
    main()
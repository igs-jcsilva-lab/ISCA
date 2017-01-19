

# This script follows global_alignment.py and shares some of the same inputs.
# The output will be a stats file describing things like the proportion of
# assembled reads that best mapped to which reference isolates. 
#
# Run the script using a command like this:
# python3 assess_alignment.py -ffs_map /path/to/format_for_spades.tsv -ga_stdout ga_out.tsv -algn_path -/path/to/alignments_out -out /path/to/stats.txt
#
# Author: James Matsumura

import re,argparse,os
from Bio import AlignIO

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    parser.add_argument('-ffs_map', type=str, required=True, help='Path to map.tsv output from format_for_spades.py.')
    parser.add_argument('-ga_stdout', type=str, required=True, help='Path to where the STDOUT of global_alignment.py went.')
    parser.add_argument('-algn_path', type=str, required=True, help='Path to the the directory preceding all the alignment directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('-out', type=str, required=True, help='Path to output directory for these stats.')
    args = parser.parse_args()

    # dict for counting which isolates have the best alignments
    isolate_counts = {} 
    
    # list for mapping percent ID match + potential coverage for
    # a given read, used to generate an outfile for plotting.
    id_v_cov = []

    # dict for counting how strong the %ID is for the best alignment
    percent_id = {"0<=x<10":0,"10<=x<20":0,"20<=x<30":0,"30<=x<40":0,
                "40<=x<50":0,"50<=x<60":0,"60<=x<70":0,"70<=x<80":0,
                "80<=x<90":0,"90<=x<100":0,"x=100":0}  

    # A set that will specify which directories of alignments to skip over.
    unassembled = set()

    # First identify which assemblies could not align. This is captured by the 
    # STDOUT of global_alignment.py
    if args.ga_stdout:
        with open(args.ga_stdout,'r') as align_map:
            for line in align_map:
                line = line.rstrip()
                ele = line.split('\t')
                unassembled.add(ele[0])

    # Need to iterate over the map generated from SPAdes step.
    with open(args.ffs_map,'r') as loc_map:
        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[0]

            # Leave early if we know this was unable to be assembled.
            if locus in unassembled:
                continue

            algn_dir = "{0}/{1}".format(args.algn_path,locus)
            isos,scores,ids,files,cov = ([] for i in range(5)) # reinitialize for every locus

            # If the minimum threshold is set high enough, it is possible for
            # no alignments to have been performed. Print to STDOUT in case
            # this does happen. 
            aligned = False

            # Found the alignment directory for this locus, now iterate over 
            # the final alignments and pull the best score.
            for file in os.listdir(algn_dir):
                a,b = ("" for i in range(2)) # store lengths of the trimmed alignments
                if file.endswith(".trimmed_align.txt"):

                    aligned = True

                    isolate = file.split('.')[0] # grab the isolate name
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
                print("The locus {0} could assemble but none of the contigs passed the minimum threshold chosen when running global_alignment.py".format(locus))

            best = scores.index(max(scores))
            best_iso = isos[best]
            best_id = ids[best]
            best_cov = cov[best]
            best_file = files[best]

            # extract the %ID match + coverage for plotting
            id_v_cov.append("{0}:{1}:{2}".format(best_id,best_cov,best_file))

            # Now assess overall stats for which reference aligned and bin the % ID. 
            if best_iso in isolate_counts:
                isolate_counts[best_iso] += 1
            else:
                isolate_counts[best_iso] = 1 
            percent_id = bin_percent_id(percent_id,best_id)

    # Write out the overall alignment stats here. 
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
        for bin,count in percent_id.items():
            rel = float("{0:.2f}".format(100*count/tot))
            if bin != 'x=100' and bin != '0<=x<10':
                out.write("{0}\t\t{1} ({2}%)\n".format(bin,count,rel))
            else:
                out.write("{0}\t\t\t{1} ({2}%)\n".format(bin,count,rel))

    # Write out a file that can generate a plot of %ID v coverage
    outfile = "{0}/ids_v_cov.tsv".format(args.out)
    with open(outfile,'w') as out:
        for pair in id_v_cov:
            ele = pair.split(':')
            out.write("{0}\t{1}\t{2}\n".format(ele[0],ele[1],ele[2]))


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
            elif line.startswith('a.trimmed'): # reached actual alignment, no need
                break

    return stats

# Function to house what is essentially a switch statement for grouping together
# different percent identity values for output.
# Argument:
# id_dict = a percent_id dict that will be built upon
# percent = percent value to use to bin
def bin_percent_id(id_dict,percent):

    percent = int(float(percent)) # can trim decimals due to bin

    # Take advantage of the order of processing to slim down the if statements
    if percent == 100:
        id_dict['x=100'] += 1
    elif percent >= 90:
        id_dict['90<=x<100'] += 1
    elif percent >= 80:
        id_dict['80<=x<90'] += 1
    elif percent >= 70:
        id_dict['70<=x<80'] += 1
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



# A script that follows threaded_assess_alignment.py and outputs some general
# statistics like counts for which reference had the best alignment per locus
# and the distribution of the %ID from the final alignments. The input for
# this script is the ids_v_cov.tsv file generated from the assessment script.
#
# Run the script using a command like this:
# python3 generate_alignment_stats.py -i ids_v_cov.tsv -o stats.out -f_or_s F
#
# Author: James Matsumura

import argparse,collections

def main():

    parser = argparse.ArgumentParser(description='Script to generate basic stats from the output of threaded_assess_alignment.py.')
    parser.add_argument('-i', type=str, required=True, help='Path to ids_v_cov.tsv output from threaded_assess_alignment.py.')
    parser.add_argument('-o', type=str, required=True, help='Name of an outfile.')
    parser.add_argument('-F_or_S', type=str, required=True, help='F or S for First or Second alignment assessment output.')
    args = parser.parse_args()

    # dict for counting which isolates have the best alignments
    isolate_counts = {} 

    # dict for counting how strong the %ID is for the best alignment
    percent_id = {"0<=x<10":0,"10<=x<20":0,"20<=x<30":0,"30<=x<40":0,
                "40<=x<50":0,"50<=x<60":0,"60<=x<70":0,"70<=x<80":0,
                "80<=x<90":0,"90<=x<100":0,"x=100":0}  
    
    with open(args.i,'r') as infile:
        for line in infile:
            line = line.rstrip()
            elements = line.split('\t')

            # Sort the %ID into bins
            id = ""
            if args.F_or_S == "F":
                id = elements[0]
            elif args.F_or_S == "S":
                id = elements[4]
            percent_id = bin_percent_id(percent_id,id)

            # Count how many times a ref is found to have the best alignment
            ref = elements[2]
            if ref in isolate_counts:
                isolate_counts[ref] += 1
            else:
                isolate_counts[ref] = 1

    # Write out the overall alignment stats here. 
    with open(args.o,'w') as out:

        out.write("Number of times a given isolate was best aligned to (if an isolate doesn't appear, never was the best alignment):\n")
        tot = 0
        for k,v in isolate_counts.items():
            tot += v # calculate a total to get relative percentage
        for k,v in isolate_counts.items():
            rel = float("{0:.2f}".format(100*v/tot)) 
            out.write("{0}\t\t{1} ({2}%)\n".format(k,v,rel))

        out.write("\nNumber of times a given percent identity was found for the best alignment:\n")
        tot = 0
        for k,v in percent_id.items():
            tot += v

        sorted_bins = collections.OrderedDict(sorted(percent_id.items(),reverse=True))

        for bin in sorted_bins:
            rel = ""
            if percent_id[bin] != 0:
                rel = float("{0:.2f}".format(100*percent_id[bin]/tot))
            else:
                rel = "0.0"
            if bin != 'x=100' and bin != '0<=x<10':
                out.write("{0}\t\t{1} ({2}%)\n".format(bin,percent_id[bin],rel))
            else:
                out.write("{0}\t\t\t{1} ({2}%)\n".format(bin,percent_id[bin],rel))


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
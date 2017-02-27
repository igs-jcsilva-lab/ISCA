

# A script that follows threaded_assess_alignment.py and outputs some general
# statistics like counts for which reference had the best alignment per locus
# and the distribution of the %ID from the final alignments. The input for
# this script is the ids_v_cov.tsv file generated from the assessment script.
#
# Run the script using a command like this:
# python3 generate_alignment_stats.py -i ids_v_cov.tsv -o stats.out
#
# Author: James Matsumura

import argparse
import numpy as np

def main():

    parser = argparse.ArgumentParser(description='Script to generate basic stats from the output of threaded_assess_alignment.py.')
    parser.add_argument('-i', type=str, required=True, help='Path to TSV output file from analyze_sequences.py/split_analysis.py')
    args = parser.parse_args()

    gc,length,rpts,rpts_len,exns = ([] for i in range(5))

    with open(args.i,'r') as infile:
        for line in infile:
            line = line.rstrip()
            elements = line.split('\t')
            gc.append(float(elements[1]))
            length.append(float(elements[2]))
            rpts.append(float(elements[3]))
            rpts_len.append(float(elements[4]))
            exns.append(float(elements[5]))

    print("GC avg:\t{0}".format(np.mean(gc)))
    print("Length avg:\t{0}".format(np.mean(length)))
    print("Repeats avg:\t{0}".format(np.mean(rpts)))
    print("Repeats length avg:\t{0}".format(np.mean(rpts_len)))
    print("Exons avg:\t{0}".format(np.mean(exns)))

    print("GC std:\t{0}".format(np.std(gc)))
    print("Length std:\t{0}".format(np.std(length)))
    print("Repeats std:\t{0}".format(np.std(rpts)))
    print("Repeats length std:\t{0}".format(np.std(rpts_len)))
    print("Exons std:\t{0}".format(np.std(exns))) 


if __name__ == '__main__':
    main()
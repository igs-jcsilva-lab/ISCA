

# This script is supplementary to the pipeline and allows one to generate
# a histogram off the TSV output file from assess_alignment.py.
#
# Run the script using a command like this:
# python3 plot_assessment.py -i aa_out.tsv
#
# Author: James Matsumura

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Script to generate a histogram of the coverage ratio for best alignments, follows assess_alignment.py.')
parser.add_argument('-i', type=str, required=True, help='Location of ids_v_coverage.tsv.')
args = parser.parse_args()

cov = []

# Just extract the coverage ratio
with open(args.i,'r') as infile:
    for line in infile:
        line = line.rstrip()
        ele = line.split('\t')
        ratio = float(ele[1])
        if ratio > 5: # trim outliers here, notify if found
            print "Potential outlier value found: %s" % (ratio)
            continue
        else:
            cov.append(float("%.2f" % ratio))

plt.hist(cov, bins=100)
plt.title('Histogram of coverage ratio (reference / assembled)')
plt.ylabel('frequency')
plt.xlabel('(reference / assembled) length ratio')
plt.show()

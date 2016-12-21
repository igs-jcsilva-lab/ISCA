

# This script is supplementary to the pipeline and allows one to generate
# a scatter plot based off the data in one of the ids_v_gaps.tsv file 
# yielded by assess_alignment.py 
#
# Run the script using a command like this:
# python3 plot_assessment.py -i ids_v_gaps.tsv
#
# Author: James Matsumura

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Script to generate a scatter plot for best alignments, follows assess_alignment.py.')
parser.add_argument('-i', type=str, required=True, help='Location of ids_v_gaps.tsv.')
args = parser.parse_args()

x,y = ([] for i in range(2))

# Assign column 1 of file to x, and column 2 to y
with open(args.i,'r') as infile:
    for line in infile:
        line = line.rstrip()
        ele = line.split('\t')
        x.append(ele[0])
        y.append(ele[1])

plt.scatter(x, y, s=10, alpha=0.5, marker="x")
plt.title('%ID vs %gaps from best alignments')
plt.ylabel('% ID match')
plt.xlabel('(reference / assembled) length ratio')
plt.ylim([0,101])
plt.yticks(range(0,101,5))
plt.show()

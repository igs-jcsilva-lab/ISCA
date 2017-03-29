

# This script is supplementary to the pipeline and allows one to generate
# a scatterplot that compares the percent identity obtained by global alignment
# versus the percent ID of the reference bases correctly accounted for. This
# takes the output from second_threaded_assess_alignment.py.
#
# Run the script using a command like this:
# python generate_repeat_plot.py -i aa_out.tsv
#
# Author: James Matsumura

import argparse
import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Script to generate a scatter plot comparing global alignment based ID and reference based ID.')
parser.add_argument('-i', type=str, required=True, help='Location of second_ids_v_coverage.tsv.')
args = parser.parse_args()

id,rid = ([] for i in range(2))

# Just extract the coverage ratio
with open(args.i,'r') as infile:
    for line in infile:
        line = line.rstrip()
        ele = line.split('\t')
        id.append(float(ele[0]))
        rid.append(float(ele[-1]))

df = pd.DataFrame({"Alignment ID":id,"Reference ID":rid},columns=["Alignment ID","Reference ID"])

splot = sns.lmplot(x="Alignment ID",y="Reference ID",
    data=df,
    fit_reg=False,
    scatter_kws={"s":10,"marker":"D"})

axes = splot.axes
axes[0,0].set_xlim(0,105)
axes[0,0].set_ylim(0,105)
splot.axes[0][0].plot((0,100), (0,100), 'k--')
#splot.show()

splot.savefig("method2_assessment.png")

# A script that follows the end of an iteration of the pipeline and accepts
# the total ids_v_cov.tsv files (from both methods of assembly). It also takes
# in a list of particular loci to assess. This will result in a TSV file with
# 3 columns. The first is the locus ID, the second is the best alignment %ID
# to the reference, and the third is whether it was best constructed from 
# assembly method 1, assembly method 2, or failed to assemble in the pipeline. 
# An optional 4th column may be present if assembly method 2 was used, this will
# note the total alignment %ID INCLUDING gaps, which is likely to be lower than
# the value in the same respective second column as that includes gaps and the
# second assembly method has a tendency to over-assemble.
#
# Run the script using a command like this:
# python3 compare_gene_list_to_assembly_results.py -l list.txt -i ids_v_cov.tsv -o stats.out
#
# Author: James Matsumura

import argparse
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to assess a list of loci and how well they were constructed by the pipeline.')
    parser.add_argument('-i', type=str, required=True, help='Path to ids_v_cov.tsv output from threaded_assess_alignment.py.')
    parser.add_argument('-l', type=str, required=True, help='List of loci to analyze their assembly results.')
    parser.add_argument('-o', type=str, required=True, help='Name of an outfile.')
    args = parser.parse_args()

    loci = set() # track all input loci
    # track which assembler did it better, second assembly set is likely to be a subset of the first assembly set
    best_assembly = defaultdict(list)

    # Grab all the loci 
    with open(args.l,'r') as loci_list:
        for line in loci_list:
            line = line.strip()
            if '.' in line:
                loci.add(line.split('.')[0])
            else:
                loci.add(line)
            

    with open(args.i,'r') as pipeline_results:
        for line in pipeline_results:
            elements = line.strip().split('\t')
            locus = elements[3].split('/')[-2]
            percent_id = ""
            if len(elements) > 4:
                percent_id = float(elements[-1])
                gap_percent_id = float(elements[0])

                if locus in best_assembly:
                    if percent_id > best_assembly[locus][0]:
                         # just clear it, this should handle crazy cases of very long id_v_cov concats
                        best_assembly[locus][:] = []
                        best_assembly[locus].append(percent_id)
                        best_assembly[locus].append(gap_percent_id)

                else:
                    best_assembly[locus].append(percent_id)
                    best_assembly[locus].append(gap_percent_id)

            else:
                percent_id = float(elements[0])

                if locus in best_assembly:
                    if percent_id > best_assembly[locus][0]:
                        best_assembly[locus][0] = percent_id

                else:
                    best_assembly[locus].append(percent_id)

    # Write out the overall stats here. 
    with open(args.o,'w') as out:
        for locus in loci:
            if locus not in best_assembly:
                out.write("{0}\t0\tNot assembled\t0\n".format(locus))

            else:
                if len(best_assembly[locus]) > 1: # asesmbled by method 2 (HGA+SPAdes)
                    out.write("{0}\t{1}\tMethod 2\t{2}\n".format(locus,best_assembly[locus][0],best_assembly[locus][1]))
                else: # assembled by method 1 (SPAdes)
                    out.write("{0}\t{1}\tMethod 1\tNA\n".format(locus,best_assembly[locus][0]))

if __name__ == '__main__':
    main()
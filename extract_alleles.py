

# This script will parse through multiple GFF3 files and pull the alleles of each.
# It expects a file as input with each line in the file being tab-delimited 
# where the first column is the type (either reference or isolate), the second
# column is the GFF3 file, the third column is the FASTA that correlates to 
# the GFF3 file, and the fourth column is the name/prefix to designate these files to. 
# While the FASTA file is not used by this script, it is used by downstream
# scripts so formatting this file like this allows for this single file 
# to be used as input for all subsequent scripts. 
#
# The input should look like this (it MUST start with whichever you want to be the reference): 
# reference     /path/to/ref.gff3   /path/to/ref.fasta  name_of_ref
# isolate       /path/to/iso1.gff3   /path/to/iso1.fasta    name_of_iso1
# isolate       /path/to/iso2.gff3   /path/to/iso2.fasta    name_of_iso2
#
# Note there can only be one reference as this is what all other alleles will map
# to. The output will be another TSV file with the first column being the reference
# ID, the second column being source/location of this gene, and the third column
# contains the start-stop coordinates for the gene. Subsequent columns will  
# have the isolates that follow this same pattern. 
#
# The output will look like this:
# ref_id_0001   ref_loc   1-8888  iso1.ref_id_0001  iso1_loc    2-7999  iso2.ref_id_0001    iso2_loc    3-8000
# 
# Run the script using a command like this:
# python3 extract_alleles.py -i /path/to/list_input.tsv -o /path/to/outfile.tsv
#
# Author: James Matsumura

import sys,re,argparse

def main():

    parser = argparse.ArgumentParser(description='Script to map alleles across GFF3 file. Read the top of the file for more details.')
    parser.add_argument('-i', type=str, required=True, help='Path to a TSV list for references and isolates.')
    parser.add_argument('-o', type=str, required=True, help='Path to where the output TSV should go.')
    args = parser.parse_args()

    i = open(args.i,'r')
    o = open(args.o,'w')

    # dictionary where the key is the ID and the value is a list for ref/loc/coords 
    allele_map = {} 

    # Iterate over each reference/isolate
    for entry in i:
        entry = entry.rstrip()
        vals = entry.split('\t')
        gff3 = vals[1]
        name = vals[3]

        # Regardless of reference or isolate, all should be mapping to the same name
        # designated by the reference. vi 
        allele_map = parse_gff3(gff3,allele_map,name)

    # Iterate over the final hash of lists and print out a TSV
    for key,value in allele_map.items():
        vals = ('\t').join(value)
        line = "{0}\t{1}\n".format(key,vals)
        o.write(line)

def parse_gff3(file,allele_map,name):

    regex_for_name = r'.*Name=([a-zA-Z0-9_\.\-]+)'

    gff3 = open(file,'r') 

    for line in gff3:
        if line.startswith('##FASTA'): # don't care about sequences
            break
        elif line.startswith('#'): # don't care about comments or header data
            pass
        else: # within the GFF3 9-column section
            ele = line.split('\t')
            if ele[2] == 'gene': # only process if it is a gene
                source = ele[0]
                start = ele[3]
                stop = ele[4]
                strand = ele[6]
                attr_name = re.search(regex_for_name,ele[8]).group(1) # extract the name from attr that links via GMAP
                id = "{0}.{1}".format(name,attr_name) # need to make a unique ID for each group/isolate that ties back to attr name
                combined_vals = "{0}|{1}|{2}|{3}|{4}".format(source,start,stop,strand,id)
                
                if attr_name not in allele_map: # initialize if not seen before
                    allele_map[attr_name] = []

                allele_map[attr_name].append(combined_vals)

    return allele_map


if __name__ == '__main__':
    main()
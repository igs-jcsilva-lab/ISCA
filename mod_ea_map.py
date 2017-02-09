

# Run the script using a command like this:
# python3 mod_ea_map.py -original ea_map.tsv -new mod_ea_map.tsv -dupes remove_duplicates.output
#
# Author: James Matsumura

import argparse,re

def main():

    parser = argparse.ArgumentParser(description='Script to modify the extract alleles map to get rid of duplicate entries.')
    parser.add_argument('-original', type=str, required=True, help='Path to extract alleles map from extract_alleles.py.')
    parser.add_argument('-new', type=str, required=True, help='Path to where the new extract alleles map.')
    parser.add_argument('-dupes', type=str, required=True, help='Conflict output from remove_duplicates.py.')
    args = parser.parse_args()

    duplicates = set()

    # First identify all the IDs for alleles that were duplicates
    with open(args.dupes,'r') as dupe_file:
        for entry in dupe_file:
            entry = entry.rstrip()
            duplicates.add(entry)

    # Now iterate over the original extract alleles map and create a new one. 
    with open(args.new,'w') as new:
        with open(args.original,'r') as original:
            for line in original:
                line = line.rstrip()
                elements = line.split('\t')
                new_line = []

                # Grab the locus
                new_line.append(elements[0])
                
                # Now dealing with the cases of particular alleles of that
                # locus potentially being considered duplicates of another
                # sequence. 
                for j in range(1,len(elements)):
                    allele = elements[j].split('|')[4]
                    
                    # If we don't find this allele in the set of duplicates,
                    # then append the entire entry (origin,coords,strand,id)
                    if allele not in duplicates:
                        new_line.append(elements[j])

                # Need to make sure that a locus is now not devoid of any
                # alleles due to them all being removed via duplicates. 
                if len(new_line) > 1:
                    outline = ("\t").join(new_line)
                    new.write("{0}\n".format(outline))


if __name__ == '__main__':
    main()

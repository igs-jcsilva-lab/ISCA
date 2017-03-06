

# This script will parse through a particular GFF3 file given a set of FASTA
# sequences and output some basic information provided by the GFF3 like the 
# description of the gene as well as the length. 
# 
# Run the script using a command like this:
# python3 extract_unaligned_info.py -fasta /path/to/in.fsa -gff3 /path/to/file.gff3 -priority 3D7 -out /path/to/outfile.tsv
#
# Author: James Matsumura

import re,argparse,urllib

def main():

    parser = argparse.ArgumentParser(description='Script to map alleles across GFF3 file. Read the top of the file for more details.')
    parser.add_argument('-fasta', type=str, required=True, help='Path to a TSV list for references and isolates.')
    parser.add_argument('-gff3', type=str, required=True, help='.')
    parser.add_argument('-priority', type=str, required=True, help='Prefix priority.')
    parser.add_argument('-outfile', type=str, required=True, help='Output file.')
    args = parser.parse_args()

    # dictionary where the key is the ID and the value is a list for ref/loc/coords 
    references = set()
    output_list = []

    # Iterate over each reference/isolate
    with open(args.fasta,'r') as fasta:
        for line in fasta:
            if line.startswith('>'+args.priority): # only grab references that we care about
                references.add(line.split('.',1)[1].strip())

    regex_for_name = r'.*Name=([a-zA-Z0-9_\.\-]+)'
    regex_for_description = r'description=(.*)'

    with open(args.gff3,'r') as gff3:
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

                    # Description is messy as it has some odd encoding
                    isolate_description = ele[8].split(';')
                    description,aliases = ("" for i in range(2))
                    for attr in isolate_description:
                        if attr.startswith('description'):
                            description = re.search(regex_for_description,attr).group(1)
                        elif attr.startswith('Alias'):
                            aliases = attr.split('=')[1].strip()

                    id = re.search(regex_for_name,ele[8]).group(1) # extract the name from attr that links via GMAP

                    if id in references:
                        output_list.append("{0}\t{1}\t{2}\n".format(id,aliases,urllib.unquote_plus(description)))

    with open(args.outfile,'w') as out:
        for ele in output_list:
            out.write(ele)


if __name__ == '__main__':
    main()
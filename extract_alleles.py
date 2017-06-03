

# This script will parse through multiple GFF3 files and pull the alleles of each.
# It expects a file as input with each line in the file being tab-delimited 
# where the first column is the type (either reference or isolate), the second
# column is the GFF3 file, the third column is the FASTA that correlates to 
# the GFF3 file, and the fourth column is the name/prefix to designate these files to. 
# While the FASTA file is not used by this script, it is used by downstream
# scripts so formatting this file like this allows for this single file 
# to be used as input for all subsequent scripts. 
#
# The input should look like this (MUST start with whichever you want to be the reference): 
# reference     /path/to/ref.gff3   /path/to/ref.fasta  name_of_ref
# isolate       /path/to/iso1.gff3   /path/to/iso1.fasta    name_of_iso1
# isolate       /path/to/iso2.gff3   /path/to/iso2.fasta    name_of_iso2
#
# *** It is VERY important that the name_of_* column does not have periods ('.') ***
# *** This is to guarantee correct mapping later on in the pipeline. ***
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
# python3 extract_alleles.py -list /path/to/list_input.tsv -insert 500 -out /path/to/outfile.tsv
#
# Author: James Matsumura

import re,argparse
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to map alleles across GFF3 file. Read the top of the file for more details.')
    parser.add_argument('-list', type=str, required=True, help='Path to a TSV list for references and isolates.')
    parser.add_argument('-insert', type=int, required=True, help='Insert size from SRA for the reads that will be used as input.')
    parser.add_argument('-out_dir', type=str, required=True, help='Directory for where the output should go.')
    args = parser.parse_args()

    # dictionary where the key is the ID and the value is a list for ref/loc/coords 
    allele_map = {} 

    # Iterate over each reference/isolate
    with open(args.list,'r') as i:
        for entry in i:
            entry = entry.rstrip()
            vals = entry.split('\t')
            type = vals[0]
            gff3 = vals[1]
            name = vals[3]

            # Regardless of reference or isolate, all should be mapping to the same name
            # designated by the reference. vi 
            allele_map = parse_gff3(gff3,allele_map,type,name,args.insert,args.out_dir,)

    # Iterate over the final hash of lists and print out a TSV
    out = "{0}/ea_map.tsv".format(args.out_dir)
    with open(out,'w') as o:
        for key,value in allele_map.items():
            vals = ('\t').join(value)
            line = "{0}\t{1}\n".format(key,vals)
            o.write(line)


# Arguments:
# file = GFF3 file
# allele_map = a dictionary with the reference ID/Name as the key and the values an allele tied to it
# ref_or_iso = is this a reference or an isolate? This will potentially change the ID
# name = prefix/name of isolate
# insert = size of the insert to check for overlaps
# out_dir = prefix for the output directory to write to
def parse_gff3(file,allele_map,ref_or_iso,name,insert,out_dir):

    regex_for_name = r'.*Name=([a-zA-Z0-9_\.\-]+)'
    regex_for_gmap_name = r'.*ID=([a-zA-Z0-9_\.\-]+)'

    # Build a dictionary of all the loci and their positions to look for any
    # instances of overlap which may assist in deciding whether to filter
    # or not when it comes to assigning individual reads per locus. 
    overlap_dict = defaultdict(list)
    intron_check = {}
    max_intron_length = {} # keep track, per allele, of the maximum intron length
    attr_name = ""

    with open(file,'r') as gff3:
        for line in gff3:
            if line.startswith('##FASTA'): # don't care about sequences
                if len(intron_check[attr_name]['list']) > 1 and ref_or_iso == "reference": # one last check for last gene
                    max_intron_length[attr_name]['max_list'] = calculate_max_intron_length(intron_check,attr_name)
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

                    if attr_name in intron_check and ref_or_iso == "reference": # make sure it's been initialized
                        if len(intron_check[attr_name]['list']) > 1: # only need to process if more than one exon
                            max_intron_length[attr_name]['max_list'] = calculate_max_intron_length(intron_check,attr_name)

                    attr_name = re.search(regex_for_name,ele[8]).group(1) # extract the name from attr that links via GMAP

                    id = "{0}.{1}".format(name,re.search(regex_for_gmap_name,ele[8]).group(1))

                    combined_vals = "{0}|{1}|{2}|{3}|{4}".format(source,start,stop,strand,id)
                    
                    if attr_name not in allele_map: # initialize if not seen before
                        allele_map[attr_name] = []

                    allele_map[attr_name].append(combined_vals)

                    overlap_dict[source].append("{0}:{1}:{2}".format(start,stop,id))

                    intron_check[attr_name] = {'list':[],'strand':""}
                    intron_check[attr_name]['strand'] = strand
                    max_intron_length[attr_name] = {'max_list':[],'start_pos':0}
                    max_intron_length[attr_name]['start_pos'] = ele[3]

                elif ele[2] == 'exon':
                    intron_check[attr_name]['list'].append("{0}:{1}".format(ele[3],ele[4]))
                        
    # Only do overlap and intron checks for the reference as we can't really
    # trust GMAP to capture these properly with the way it maps fragments of
    # genes. 
    if ref_or_iso == "reference":
        # Identify whether there is any overlap. If there is, print to STDOUT. 
        # Despite the nested loops, this shouldn't be too bad since it's split
        # up by each GFF3 file. Doing it this way since the GFF3 files aren't 
        # guaranteed to be in order and genes are not always going to overlap
        # in a consistent manner (some may span multiple, just 1 bp, etc.)
        overlap_set = set() # don't add duplicates
        overlap_out = "{0}/overlap.tsv".format(out_dir)
        with open(overlap_out,'a') as out:
            for k,v in overlap_dict.items():
                for j in range(0,len(v)):

                    ele = v[j].split(':')
                    jstart = int(ele[0])
                    jstop = int(ele[1])
                    jid = ele[2]

                    for x in range(0,len(v)):
                        if j != x: # only compare to other gene regions
                            ele = v[x].split(":")
                            xstart = int(ele[0])
                            xstop = int(ele[1])
                            xid = ele[2]
                            pair = ""

                            if xid < jid:
                                pair = "{0}{1}".format(xid,jid)
                            else:
                                pair = "{0}{1}".format(jid,xid)
                                
                            if jstart < (xstart-insert) < jstop:
                                if pair not in overlap_set:
                                    overlap_set.add(pair)
                                    out.write('{0}\t{1}\n'.format(jid,xid))
                            elif jstart < (xstop+insert) < jstop: 
                                if pair not in overlap_set:
                                    overlap_set.add(pair)
                                    out.write('{0}\t{1}\n'.format(jid,xid))

        # Write out some output for intron length
        introns_out = "{0}/intron_positions.tsv".format(out_dir)
        with open(introns_out,'a') as out:
            for k in max_intron_length:
                if max_intron_length[k]['max_list']:
                    out.write("{0}\t{1}\t{2}\n".format(k,max_intron_length[k]['start_pos'],"\t".join(max_intron_length[k]['max_list'])))

    return allele_map


def calculate_max_intron_length(intron_dict,key):

    prev = -1 # previous end position
    max = 0 # maximum intron length found so far
    out_list = [0]

    if intron_dict[key]['strand'] == '-':

        for exon in reversed(intron_dict[key]['list']):
            out_list.append(exon)
            if prev == -1:
                prev = exon.split(":")[1]
            else:
                intron_length = int(exon.split(":")[0])-int(prev)
                if intron_length > max:
                    max = intron_length
                prev = exon.split(":")[1]

    else:

        for exon in intron_dict[key]['list']:
            out_list.append(exon)
            if prev == -1:
                prev = exon.split(":")[1]
            else:
                intron_length = int(exon.split(":")[0])-int(prev)
                if intron_length > max:
                    max = intron_length
                prev = exon.split(":")[1]
    
    out_list[0] = str(max)           
    return out_list


if __name__ == '__main__':
    main()
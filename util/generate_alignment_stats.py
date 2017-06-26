

# A script that follows threaded_assess_alignment.py and outputs some general
# statistics like counts for which reference had the best alignment per locus
# and the distribution of the %ID from the final alignments. The input for
# this script is the ids_v_cov.tsv file generated from the assessment script.
#
# Run the script using a command like this:
# python3 generate_alignment_stats.py -i ids_v_cov.tsv -o stats.out -ea_map /path/to/extract_alleles.output -groupby l
#
# Author: James Matsumura

import argparse,collections
from collections import defaultdict
from statistics import mean

def main():

    parser = argparse.ArgumentParser(description='Script to generate basic stats from the output of threaded_assess_alignment.py.')
    parser.add_argument('-i', type=str, required=True, help='Path to ids_v_cov.tsv output from threaded_assess_alignment.py.')
    parser.add_argument('-groupby', type=str, required=True, help='Get best %s by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds".')
    parser.add_argument('-ea_map', type=str, required=True, help='Path to map.tsv output from extract_alleles.py.')
    parser.add_argument('-o', type=str, required=True, help='Name of an outfile.')
    args = parser.parse_args()

    best_id = {} # dict for capturing best ID per locus/exon
    cds_map = defaultdict(list)
    cds_lengths = {} # count how many exons in a CDS from ea_map

    # dict for counting how strong the %ID is for the best alignment
    percent_id = {"0<=x<10":0,"10<=x<20":0,"20<=x<30":0,"30<=x<40":0,
                "40<=x<50":0,"50<=x<60":0,"60<=x<70":0,"70<=x<80":0,
                "80<=x<90":0,"90<=x<100":0,"x=100":0} 

    with open(args.i,'r') as infile:
        for line in infile:
            line = line.rstrip()
            elements = line.split('\t')

            entity = ""
            if args.groupby == 'l':
                if 'exon_' in elements[3]:
                    entity = elements[3].split('/')[-2]
                else:
                    entity = elements[3].split('.')[1]
            else:
                entity = elements[3].split('/')[-1].split('.WITH')[0]

            # Sort the %ID into bins
            id = 0.0
            if len(elements) == 4:
                id = float(elements[0])
            else:
                id = float(elements[4])

            if entity in best_id:
                if id > best_id[entity]:
                    best_id[entity] = id

            else:
                best_id[entity] = id

    for k,v in best_id.items():
        percent_id = bin_percent_id(percent_id,v)

        if args.groupby == 'cds':
            parent = get_exon_parent(k)
            cds_map[parent].append(v)

    tot,entity_count,loci = (0 for i in range(3))

    with open(args.ea_map,'r') as infile:
        for line in infile:
            loci += 1
            if args.groupby != 'l':
                entity_count += (len(line.split('\t'))-1) # count number of alleles/exons

                if args.groupby == 'cds':
                    elements = line.split('\t')
                    for x in range(1,len(elements)):
                        parent = get_exon_parent(elements[x].split('|')[-1]).strip()
                        if parent in cds_lengths:
                            cds_lengths[parent] += 1
                        else:
                            cds_lengths[parent] = 1

    if args.groupby == 'l':
        entity_count = loci

    # Write out the overall alignment stats here. 
    with open(args.o,'w') as out:

        for k,v in percent_id.items():
            tot += v

        out.write("\nTotal number of entities (chosen by -groupby variable):\t{1}\n".format(entity,entity_count))
        out.write("\nNumber of times a given percent identity was found for the best alignment:\n")

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

        failures = entity_count-tot
        out.write("{0}\t\t\t{1} ({2}%)\n".format("Failed to align/assemble",failures,float("{0:.2f}".format(100*failures/tot))))

        if args.groupby == 'cds':

            out.write("\nNumber of times a CDS had the following exon percent identities (each x is an exon):\n")

            cds_bins = {"all x = 100%":0,"all 90%<=x<=100%":0,"some x >= 90%, some x< 90%":0,"all < 90%":0,"partial":0}

            for k,v in cds_map.items():

                if cds_lengths[k] != len(v):
                    cds_bins["partial"] += 1
                else:
                    cds_bins = bin_cds_stats(cds_bins,v)

            for k,v in cds_bins.items():
                out.write("\n{0}\t\t\t{1}\n".format(k,v))

            
def get_exon_parent(exon):
    if '-' in exon or 'exon_' in exon:
        return exon.split('-')[0]
    else:
        return exon.rsplit('.',1)[0]


# Function to house what is essentially a switch statement for grouping together
# different percent identity values for CDS output.
# Argument:
# cds_bins = dict for bins 
# exon_list = a list consisting of exons that comprise a single CDS 
def bin_cds_stats(cds_bins,exon_list):

    passed_threshold = True
    for exon in exon_list:
        if int(exon) == 100:
            continue
        else:
            passed_threshold = False
            break
    if passed_threshold:
        cds_bins["all x = 100%"] += 1
        return cds_bins


    passed_threshold = True
    for exon in exon_list:
        if int(exon) >= 90:
            continue
        else:
            passed_threshold = False
            break
    if passed_threshold:
        cds_bins["all 90%<=x<=100%"] += 1
        return cds_bins

    passed_threshold = True
    for exon in exon_list:
        if int(exon) < 90:
            continue
        else:
            passed_threshold = False
            break
    if passed_threshold:
        cds_bins["all < 90%"] += 1
        return cds_bins

    else:
        cds_bins["some x >= 90%, some x< 90%"] += 1
        return cds_bins

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
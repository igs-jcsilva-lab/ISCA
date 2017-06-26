

# A script that follows threaded_assess_alignment.py and pulls the best 
# assembled sequence into a single FASTA file. 
#
# Run the script using a command like this:
# python3 generate_FASTA.py -ivc ids_v_cov.tsv -outfile out.fsa -threshold 90
#
# Author: James Matsumura

import argparse,collections
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description='Script to generate basic stats from the output of threaded_assess_alignment.py.')
    parser.add_argument('-ivc', type=str, required=True, help='Path to ids_v_cov.tsv output from threaded_assess_alignment.py.')
    parser.add_argument('-threshold', type=int, default=0, required=False, help='Cutoff for pulling a sequence or not. Can set to 0 to get anything.')
    parser.add_argument('-groupby', type=str, required=True, help='Get sequences by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds".')
    parser.add_argument('-outfile', type=str, required=True, help='Name of an outfile.')
    args = parser.parse_args()

    best_id,cds_map = (defaultdict(list) for i in range(2)) 
    cds_lengths = {} # count how many exons in a CDS from ea_map


    with open(args.ivc,'r') as infile:
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

            if int(id) < args.threshold:
                continue

            if entity in best_id:
                if id > best_id[entity][0]:
                    best_id[entity][0] = id
                    best_id[entity][1] = elements[3]

            else:
                best_id[entity].append(id)
                best_id[entity].append(elements[3])

    if args.groupby == 'cds':
        for k,v in best_id.items():
            parent = get_exon_parent(k)
            cds_map[parent].append(v)

    final_sequences = []
    for k,v in best_id.items():
        new_id = "assembled_{0}".format(k)

        file = v[1].replace('trimmed_align.txt','b.fsa')
        record = SeqIO.read(file, "fasta")

        record.id = new_id
        record.description = ''
        record.name = ''
        if '.r.trimmed' in v[1]:
            tmp_seq = Seq(str(record.seq))
            record.seq = tmp_seq.reverse_complement() 

        final_sequences.append(record)

    SeqIO.write(final_sequences, args.outfile, 'fasta')
            
def get_exon_parent(exon):
    if '-' in exon or 'exon_' in exon:
        return exon.split('-')[0]
    else:
        return exon.rsplit('.',1)[0]

if __name__ == '__main__':
    main()
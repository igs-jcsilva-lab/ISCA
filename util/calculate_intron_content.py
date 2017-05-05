

# This script follows extract_sequences.py and reports the AT content of the
# intronic region for each sequence (if introns are present). It does this
# by pooling together all the introns and considers those one long sequence
# where it calculates AT content of that sequence.
#
# Run the script using a command like this:
# python3 calculate_intron_content.py -intron_map intron_map -seqs seqs.fsa -outfile output
#
# Author: James Matsumura

import re,argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

def main():

    parser = argparse.ArgumentParser(description='Script to assess the results of the base Targeted Assembly pipeline.')
    parser.add_argument('-intron_map', type=str, required=True, help='Path .')
    parser.add_argument('-seqs', type=str, required=True, help='FASTA file that corresponds to the introns in the map.')
    parser.add_argument('-min_intron_length', type=int, required=True, help='Min length of intron to keep in output.')
    parser.add_argument('-final_verdict_map', type=str, required=True, help='Map derived from the final verdict script.')
    parser.add_argument('-outfile', type=str, required=True, help='Path to where the unaligned/unassembled FASTA entries and the new alignments map should go.')
    args = parser.parse_args()

    seq_dict = SeqIO.to_dict(SeqIO.parse(args.seqs,"fasta"))
    best_id = {}

    with open(args.final_verdict_map,'r') as fvm:
        for line in fvm:
            result = line.strip().split('\t')
            regex_for_locus = r'/([a-zA-Z0-9\_\.]+).txt'
            locus = re.search(regex_for_locus,result[3]).group(1).split('.')[1]
            if locus in best_id:
                if float(result[0]) > float(best_id[locus]):
                    best_id[locus] = result[0]
            else:
                best_id[locus] = result[0]

    with open(args.intron_map,'r') as im:
        with open(args.outfile,'w') as out:
            out.write("locus\toverall_gc\texon_gc\tintron_gc\tmax_intron_length\tgene_length\tlength_after_splicing\tbest_id\n")
            for line in im:
                elements = line.strip().split("\t")

                # already re-oriented these in extract_alleles.py
                locus = "3D7.{0}".format(elements[0])
                mapped_locus = locus.split('.')[1]
                start_pos = elements[1]
                exons = elements[3:len(elements)]
                max_intron_length = int(elements[2])

                if args.min_intron_length:
                    if max_intron_length < args.min_intron_length:
                        continue

                intron_stats = calculate_intron_content(start_pos,seq_dict[locus].seq,exons)
                exon_only = calculate_exon_content(start_pos,seq_dict[locus].seq,exons)
                seq_len = len(seq_dict[locus].seq)
                exon_only_length = int(seq_len) - int(intron_stats[0])
                id = "0"
                if mapped_locus in best_id:
                    id = best_id[mapped_locus]
                
                out.write("{0}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4}\t{5}\t{6}\t{7}\n".format(elements[0],GC(seq_dict[locus].seq),exon_only,intron_stats[1],max_intron_length,seq_len,exon_only_length,id))

# Feed this the start_pos/offset for the intron positions and the sequence
# to pull introns from and it gives back GC content for just the intronic 
# sequences.
def calculate_intron_content(start_pos,seq,exons):

    intron_seq = ""

    for j in range(0,len(exons)-1):
        exon1_end = int(exons[j].split(':')[1]) - int(start_pos)
        exon2_sta = int(exons[j+1].split(':')[0]) - int(start_pos)

        intron_seq += seq[exon1_end:exon2_sta]

    intron_details = []
    intron_details.append(len(intron_seq))
    intron_details.append(GC(intron_seq))
    return intron_details

def calculate_exon_content(start_pos,seq,exons):

    exon_seq = ""

    for exon in exons:
        sta = int(exon.split(':')[0]) - int(start_pos)
        end = int(exon.split(':')[1]) - int(start_pos)

        exon_seq += seq[sta:end]

    return GC(exon_seq)

if __name__ == '__main__':
    main()
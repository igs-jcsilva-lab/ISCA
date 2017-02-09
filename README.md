# targeted_assembly
A new pipeline to perform targeted assembly given WGS reads, reference genome assemblies, and a primary reference annotation (GFF3)

# Complete steps:
1. Map an annotated reference genome to other assembled genomes
- GMAP
2. Build a map for the alleles extracted from GFF3
-  extract_alleles.py 
3. Extract sequences for all references given the previous scripts output
- extract_sequences.py 
4. remove duplicate FASTA sequences
- remove_duplicates.py
4.1 If duplicates are found, need to reformat the map from extract_alleles.py to not contain these sequences.
- mod_ea_map.py
5. Bowtie2
5.1 Build index
5.2 align
5.3 optional, but can compress SAM to BAM here
6. Analyze BAM to map reads to refs and vice-versa 
- analyze_bam.py
7. Assign all the reads to their own directories for each reference
- fastq_reads_to_fastq_alleles.py
8. Rename all the directories to format for running SPAdes on the grid 
- format_for_spades.py
9. SPAdes
10. Run global alignment 
- threaded_global_alignment.py
11. Run assessment to isolate the best assemblies and overall stats
- assess_alignment.py

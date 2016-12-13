# targeted_assembly
A new pipeline to perform targeted assembly given WGS reads, reference genome assemblies, and a primary reference annotation (GFF3)

# Order of the scripts:
1. extract_alleles.py
2. extract_sequences.py
3. analyze_bam.py
4. fastq_reads_to_fastq_alleles.py
5. format_for_spades.py

# Complete steps:
1. GMAP
2. extract_alleles.py
3. extract_sequences.py
4. remove duplicate FASTA sequences
5. Bowtie2
5.1 Build index
5.2 align
6. analyze_bam.py
7. fastq_reads_to_fastq_alleles.py
8. format_for_spades.py
9. SPAdes
10. EMBOSS

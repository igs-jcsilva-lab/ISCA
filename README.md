# targeted_assembly
A pipeline to perform targeted assembly of individual loci given WGS reads, reference genome assemblies, and a primary reference annotation (GFF3)

# Complete steps:
1. Map an annotated reference genome to other assembled genomes
  * GMAP
2. Build a map for the alleles extracted from GFF3
  * extract_alleles.py 
3. Extract sequences for all references given the previous scripts output
  * extract_sequences.py 
4. remove duplicate FASTA sequences
  * remove_duplicates.py
  1. If duplicates are found, need to reformat the map from extract_alleles.py to not contain these sequences.
    * mod_ea_map.py
5. Bowtie2
  1. Build index
  2. align
  3. optional, but can compress SAM to BAM here
6. Analyze BAM to map reads to refs and vice-versa 
  * analyze_bam.py
7. Assign all the reads to their own directories for each reference
  * fastq_reads_to_fastq_alleles.py
8. Rename all the directories to format for running SPAdes on the grid 
  * format_for_spades.py
9. SPAdes
  * http://spades.bioinf.spbau.ru/release3.5.0/manual.html
10. Run global alignment 
  * threaded_global_alignment.py
11. Run assessment to isolate the best assemblies and overall stats
  * threaded_assess_alignment.py
12. If there are any remaining loci that could not assemble at a desired minimum threshold, can isolate these reference sequences to another round of the pipeline and use more sensitive Bowtie2 alignment parameters. Note that using this step will essentially format the data similar to the end of step 3. 
  * extact_new_round_seqs.py

Within the util/~ directory there are a number of other post-processing scripts that can
be used to analyze the results of this final step. 
* generate_alignment_stats.py - will give an overview of the output of threaded_assess_alignment.py
* generate_histogram.py - will generate a histogram plot of the coverage found from the output of threaded_assess_alignment.py
* generate_scatterplot.py - will generate a scatter plot (%ID v coverage) from the output of threaded_assess_alignment.py


Dependencies:
- Python 3
  * [Biopython](https://pypi.python.org/pypi/biopython/1.66)
    * [EMBOSS](http://emboss.open-bio.org/)
  * [pysam](https://pypi.python.org/pypi/pysam)
- [SPAdes 3.5+](http://bioinf.spbau.ru/spades)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

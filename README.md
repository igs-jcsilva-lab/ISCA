# targeted_assembly
A pipeline to perform targeted assembly of individual loci given WGS reads, reference genome assemblies, and a primary reference annotation (GFF3)

## Dependencies
- Python 3.5
  * [Biopython](https://pypi.python.org/pypi/biopython/1.66)
  * [pysam](https://pypi.python.org/pypi/pysam)
- Python 2.7 (Needed for the externally developed scripts, HGA+Scaffold Builder, as well as CWL)
  * [cwlref-runner](https://pypi.python.org/pypi/cwlref-runner)
  * [pyyaml](https://pypi.python.org/pypi/PyYAML)
- [GSNAP](http://research-pub.gene.com/gmap/) - tested with release 2017-01-14
- [SMALT](http://www.sanger.ac.uk/science/tools/smalt-0) - tested with v0.7.6
- [SPAdes](http://bioinf.spbau.ru/spades) - tested with v3.10.1
- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) - tested with v1.2.10
- [EMBOSS](http://emboss.open-bio.org/) - tested with v6.6.0.0
- Python Scripts (automatically pulled when running the pipeline)
  * [Hierarchical Genome Assembly Tool (HGA)](https://github.com/jmatsumura/Hierarchical-Genome-Assembly-HGA)
    * [Original for reference](https://github.com/aalokaily/Hierarchical-Genome-Assembly-HGA)
  * [Scaffold Builder](https://github.com/jmatsumura/Scaffold_builder)
    * [Original for reference](https://github.com/metageni/Scaffold_builder)

## Required Inputs
* Annotated GFF3+FASTA file for a genome (FASTA containing the whole chromosomes/sequence-regions that matches the positions noted in GFF3)
* OPTIONAL 
  * A mapped (via GMAP) annotated reference genome to any others if one wants to pool together reference loci to recruit reads for one locus' assembly across alleles
* A TSV file noting the path of the GFF3+FASTA file(s) in the format specified in the top of [extract_alleles.py](https://github.com/jmatsumura/targeted_assembly/blob/master/bin/extract_alleles.py) -- (example)[https://github.com/jmatsumura/targeted_assembly/blob/master/example_data/ea_input.tsv]
* A list of loci to focus the assembly on -- (example)[https://github.com/jmatsumura/targeted_assembly/blob/master/example_data/subset_list.txt]
* Two FASTQ paired reads files from the isolate to perform targeted assembly with

## Pipeline steps (all scripts found in ./bin)
1. Build a map for the alleles extracted from GFF3
  * `extract_alleles.py` 
2. Extract sequences for all references given the previous scripts output
  * `extract_sequences.py`
3. Set up the working directories and locations for outputs
  * `build_workspace.py`
4. GSNAP
  * Build index
  * align
  * optional, but can compress SAM to BAM here
5. Analyze BAM to map reads to refs and vice-versa 
  * `analyze_bam.py`
6. Assign all the reads to their own directories for each reference
  * `fastq_reads_to_fastq_alleles.py`
7. Rename all the directories to format for running SPAdes on the grid 
  * `format_for_assembly.py`
8. Assembly (either via SPAdes or HGA)
  * `run_parallel_assembly.py`
9. Run alignment 
  * `threaded_alignment.py`
10. Run assessment to isolate the best assemblies and overall stats
  * `threaded_assess_alignment.py`
11. If there are any remaining loci that could not assemble at a desired minimum threshold, can isolate these reference sequences to another round of the pipeline and use a different aligner/sensitivity. Note that using this step will essentially format the data similar to the end of step 3
  * `assembly_verdict.py`
12. Try assemble those that the previous assembler could not
  * `run_parallel_assembly.py`
13. Rerun alignment using these new assemblies.
  * `threaded_alignment.py`
14. Assess these new assemblies.
  * `threaded_assess_alignment.py`
15. Build a dataset for those that cannot align
  * `assembly_verdict.py`
16. Repeat steps 4-15, but at step 4 use the SMALT aligner

## Invoking a CWL workflow
```
cwl-runner <cwl tool/workflow script> <input parameter yml/json>
```
The first parameter is a valid cwl tool or workflow script.  These have the extension __.cwl__.

The second parameter is a YAML or JSON file consisting of input parameters for the CWL script. YAML examples are provided and are listed with the extension __.yml__.
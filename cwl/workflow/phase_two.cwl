#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Phase two of the workflow consists of analyze_bam.cwl, fastq_reads_to_fastq_alleles.cwl, and format_for_assembly.cwl
class: Workflow


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  sam:
    label: Path to a SAM file derived from aligner
    type: File

  ea_map:
    label: Path to the output from extract_alleles.py
    type: File

  threshold:
    label: Percent cut-off for how many matches a read must hit (uses whichever is longer the read or reference as denominator)
    type: int

  prefix:
    label: Name of the prefix to yield the two maps (one read-based and one reference-based)
    type: string

  reads_dir:
    label: Reads directory made by build_workspace.cwl
    type: Directory

  filter:
    label: Either "yes" or "no" for removing discrepancies + multi-locus mapping reads
    type: string

  paired_suffixes:
    label: Either "yes" or "no" for whether the reads are mapped to one another with suffixes like .1 and .2 and one wants to assess for concordancy. This is dependent on the aligner. Check the *read_map.tsv file and see if the first elements are by read pair (so no suffix) or individual read (each read has  suffix) and answer accordingly
    type: string

  reads1:
    label: Path to where the first paired fastq.gz file
    type: File

  reads2:
    label: Path to where the second paired fastq.gz file
    type: File

  assmb_path:
    label: Path to the the directory to initialize directories for all the assembly output
    type: Directory

  outfile:
    label: Name of the map to create which maps a reference locus to an int ID
    type: string

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  read_map:
    type: File
    outputSource: analyze_bam/read_map

  ref_map:
    type: File
    outputSource: analyze_bam/ref_map

  renamed_reads_dir:
    type: Directory
    outputSource: format_for_assembly/renamed_reads_dir

  renamed_assmb_dir:
    type: Directory
    outputSource: format_for_assembly/renamed_assmb_dir

  assmb_map:
    type: File
    outputSource: format_for_assembly/assmb_map


steps:
  analyze_bam:
    run: analyze_bam.cwl
    in:
      sam: sam
      ea_map: ea_map
      threshold: threshold
      prefix: prefix
      python3_lib: python3_lib
    out: [read_map,ref_map]

  fastq_reads_to_fastq_alleles:
    run: fastq_reads_to_fastq_alleles.cwl
    in:
      reads_dir: reads_dir
      filter: filter
      paired_suffixes: paired_suffixes
      ab_read_map: analyze_bam/read_map
      reads1: reads1
      reads2: reads2
      python3_lib: python3_lib
    out: [assigned_reads_dir]
    
  format_for_assembly:
    run: format_for_assembly.cwl
    in:
      reads_dir: fastq_reads_to_fastq_alleles/assigned_reads_dir
      assmb_path: assmb_path
      ab_ref_map: analyze_bam/ref_map
      outfile: outfile
      python3_lib: python3_lib
    out: [renamed_reads_dir,renamed_assmb_dir,assmb_map]
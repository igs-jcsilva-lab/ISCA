#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: GSNAP - complete workflow for building a GSNAP index and then aligning the reads to the index
class: Workflow


requirements:
  - class: InlineJavascriptRequirement

inputs:
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string

  gsnap_dir:
    label:  Path to the output directory to write the GSNAP genome files to
    type: Directory

  sequences:
    label: Path to the sequence file built from extract_sequences.py
    type: File

  threads:
    label: Number of threads to use for alignment
    type: int

  reads1:
    label: Path to the first read pair file
    type: File

  reads2:
    label: Path to the second read pair file
    type: File

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  gsnap_sam:
    type: File
    outputSource: gsnap_align/gsnap_sam


steps:
  gsnap_index:
    run: gsnap_index.cwl
    in:
      gsnap_genome: gsnap_genome
      gsnap_dir: gsnap_dir
      sequences: sequences
      python3_lib: python3_lib
    out: [stdout]

  gsnap_align:
    run: gsnap_align.cwl
    in:
      gsnap_genome: gsnap_genome
      gsnap_dir: gsnap_dir
      threads: threads
      reads1: reads1
      reads2: reads2
      python3_lib: python3_lib
      gsnap_index_stdout: gsnap_index/stdout
    out: [gsnap_sam]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: GSNAP - complete workflow for building a GSNAP index and then aligning the reads to the index
class: Workflow


requirements:
  - class: InlineJavascriptRequirement


inputs:
# Shared values
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string

  gsnap_dir:
    label: Path to the output directory to write the GSNAP genome files to
    type: string

  sequences:
    label: Path to the sequence file built from extract_sequences.py
    type: string

  threads:
    label: Number of threads to use for alignment
    type: int

  reads1:
    label: Path to the first read pair file
    type: File

  reads2:
    label: Path to the second read pair file
    type: File


outputs:
  gsnap_gsnap:
    type: File
    outputSource: gsnap_align/gsnap_sam


steps:
  gsnap_index:
    run: gsnap_index.cwl
    in:
      gsnap_genome: gsnap_genome
      gsnap_dir: gsnap_dir
      sequences: sequences
    out: [gsnap_index_dir]

  gsnap_align:
    run: gsnap_align.cwl
    in:
      gsnap_genome: gsnap_genome
      gsnap_dir: gsnap_index/gsnap_index_dir
      threads: threads
      reads1: reads1
      reads2: reads2
    out: [gsnap_sam]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: SMALT - complete workflow for building a SMALT index and then aligning the reads to the index
class: Workflow


requirements:
  - class: InlineJavascriptRequirement


inputs:
  smalt_prefix:
    label: Name of the prefix used by SMALT
    type: string

  sequences:
    label: Reference sequence to build the index from
    type: File

  smalt_dir:
    label: Path to the base directory where the SMALT index files will be placed
    type: Directory

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
  smalt_files:
    type: 
      type: array
      items: File
    outputSource: smalt_index/smalt_index_files
  smalt_sam:
    type: File
    outputSource: smalt_align/smalt_sam


steps:
  smalt_index:
    run: smalt_index.cwl
    in:
      smalt_dir: smalt_dir
      smalt_prefix: smalt_prefix
      sequences: sequences
    out: [smalt_index_files]

  smalt_align:
    run: smalt_align.cwl
    in:
      smalt_dir: smalt_dir
      smalt_prefix: smalt_prefix
      threads: threads
      reads1: reads1
      reads2: reads2
    out: [smalt_sam]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: GSNAP - Build an index
class: CommandLineTool
stdout: gsnap_index.out

requirements:
  - class: InlineJavascriptRequirement


inputs:
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string
    inputBinding:
      prefix: "-d"
      position: 1

  gsnap_dir:
    label: Path to the output directory to write the GSNAP genome files to
    type: Directory
    inputBinding:
      prefix: "-D"
      position: 2

  sequences:
    label: Path to the sequence file built from extract_sequences.py
    type: File
    inputBinding:
      position: 3

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  stdout:
    type: stdout


baseCommand: ["GMAP_GSNAP_BIN/gmap_build"]
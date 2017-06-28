#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: SMALT - Build an index
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  kmer:
    label: kmer size to break the reference sequences into
    type: int
    inputBinding:
      prefix: "-k"
      position: 1
    default: 13

  step:
    label: step size for the kmers
    type: int
    inputBinding:
      prefix: "-s"
      position: 2
    default: 2

  smalt_prefix:
    label: Path to directory to build the SMALT index
    type: string
    inputBinding:
      position: 3

  sequences:
    label: Path to the sequence file built from extract_sequences.py (or assembly_verdict.py)
    type: File
    inputBinding:
      position: 4


outputs:
  smalt_index_dir:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.smalt_prefix + "*")


baseCommand: ["/usr/local/packages/smalt/bin/smalt","index"]
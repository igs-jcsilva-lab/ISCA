#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Generate sequences from the allele map
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)


inputs:
  ea_input:
    label: Path to a TSV list for references and isolates
    type: File
    inputBinding:
      prefix: "--ea_input"

  ea_map:
    label: Path to the output from extract_alleles.py
    type: File
    inputBinding:
      prefix: "--ea_map"

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int
    inputBinding:
      prefix: "--buffer"

  subset_list:
    label: List to subset the sequences by
    type: File
    inputBinding:
      prefix: "--subset_list"

  prefix:
    label: Prefix for the output FASTA files to generate in current or existing directory
    type: string
    inputBinding:
      prefix: "--prefix"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  buffered_sequences:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*_buffered*")

  unbuffered_sequences:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*_unbuffered*")


baseCommand: ["PYTHON3_BIN/python","TARGETED_ASSEMBLY_BIN/extract_sequences.py"]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Generate allele map
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


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

  outfile:
    label: Name of the output FASTA file to generate in current or existing directory
    type: string
    inputBinding:
      prefix: "--outfile"


outputs:
  sequences:
    type: File
    outputBinding:
      glob: $(inputs.outfile)


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/extract_sequences.py"]
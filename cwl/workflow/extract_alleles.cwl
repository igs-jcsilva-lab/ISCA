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

  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string
    inputBinding:
      prefix: "--gene_or_exon"


outputs:
  ea_map:
    type: File
    outputBinding:
      glob: $('ea_map.tsv')


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/extract_alleles.py"]
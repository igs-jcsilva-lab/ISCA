#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Generate allele map
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  ea_input:
    inputBinding:
      prefix: "-ea_input"
    label: Path to a TSV list for references and isolates
    type: File
  out_dir:
    inputBinding:
      prefix: "-out_dir"
    label: Location to write out the various files from extract_alleles.py
    type: string
    default: '.'
  gene_or_exon:
    inputBinding:
      prefix: "-gene_or_exon"
    label: Either "gene" or "exon" for which sequences to pull
    type: string


outputs:
  ea_map:
    type: File
    outputBinding:
      glob: $(inputs.out_dir + '/ea_map.tsv')


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/targeted_assembly/bin/extract_alleles.py"]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Generate allele map
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

  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string
    inputBinding:
      prefix: "--gene_or_exon"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  ea_map:
    type: File
    outputBinding:
      glob: $('ea_map.tsv')


baseCommand: ["PYTHON3_BIN/python","TARGETED_ASSEMBLY_BIN/extract_alleles.py"]
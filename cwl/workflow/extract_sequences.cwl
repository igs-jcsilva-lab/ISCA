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


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/extract_sequences.py"]
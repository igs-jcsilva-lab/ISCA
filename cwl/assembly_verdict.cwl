#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- assign reads to individual locus directories
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  ivc:
    label: Path to an ids_v_cov.tsv file from the previous run
    type: File
    inputBinding:
      prefix: "--ivc"

  threshold:
    label: Minimum threshold of %ID that needs to be met to pass final assembly
    type: int
    inputBinding:
      prefix: "--threshold"

  original_buffered_fsa:
    label: Path to where the unbuffered FASTA file generated from extract_sequences is
    type: File
    inputBinding:
      prefix: "--original_buffered_fsa"

  original_assmb_map:
    label: Path to where the output from format_for_assembly.py is located
    type: File
    inputBinding:
      prefix: "--original_assmb_map"

  prefix:
    label: Prefix for the output FASTA files to generate in current or existing directory
    type: string
    inputBinding:
      prefix: "--prefix"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?

  original_fsa:
    label: Path to where the unbuffered FASTA file generated from extract_sequences is
    type: File
    inputBinding:
      prefix: "--original_fsa"

  min_align_len:
    label: Minimum alignment length ratio of assembled sequence as cutoff to pull sequence or not
    type: double
    inputBinding:
      prefix: "--min_align_len"

outputs:
  leftovers:
    type: File
    outputBinding:
      glob: $("*leftovers*")

  hga_assmb_map:
    type: File
    outputBinding:
      glob: $("*assmb*")


baseCommand: ["PYTHON3_EXE","TARGETED_ASSEMBLY_BIN/assembly_verdict.py"]

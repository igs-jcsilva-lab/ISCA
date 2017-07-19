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

  original_fsa:
    label: Path to where the unbuffered FASTA file generated from extract_sequences is
    type: File
    inputBinding:
      prefix: "--original_fsa"

  original_assmb_map:
    label: Path to where the output from format_for_assembly.py is located
    type: File
    inputBinding:
      prefix: "--original_assmb_map"

  out_dir:
    label: Path to where the unaligned/unassembled FASTA entries and the new alignments map should go
    type: Directory
    inputBinding:
      prefix: "--out_dir"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  end_results:
    type: Directory
    outputBinding:
      outputEval: $(inputs.out_dir)

  hga_assmb_map:
    type: File
    outputBinding:
      glob: $("*new_assmb_map*")


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/assembly_verdict.py"]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- pull sequences reconstructed below threshold
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  workspace_location:
    label: Path where directories are built at
    type: string
    inputBinding:
      prefix: "--workspace_location"
    default: "."

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?

  first_align_path:
    label: Path to output directory for all these alignments.
    type: Directory
    inputBinding:
      prefix: "--first_align_path"

  second_align_path:
    label: Path to output directory for all these alignments.
    type: Directory
    inputBinding:
      prefix: "--second_align_path"

  first_ids_v_cov:
    label: File with metadata of assemblies that created alignments from first round
    type: File
    inputBinding:
      prefix: "--first_ids_v_cov"

  second_ids_v_cov:
    label: File with metadata of assemblies that created alignments from second round
    type: File
    inputBinding:
      prefix: "--second_ids_v_cov"

  first_idvc:
    label: File with metadata of assemblies that created alignments from first round
    type: File
    inputBinding:
      prefix: "--first_idvc"

  second_idvc:
    label: File with metadata of assemblies that created alignments from second round
    type: File
    inputBinding:
      prefix: "--second_idvc"

  original_fsa:
    label: Path to where the unbuffered FASTA file generated from extract_sequences is
    type: File
    inputBinding:
      prefix: "--original_fsa"

outputs:
  assembled_below_threshold:
    type: File
    outputBinding:
      glob: $(inputs.workspace_location + '/assembled_seqs_below_threshold.fsa')

baseCommand: ["PYTHON3_EXE","TARGETED_ASSEMBLY_BIN/extract_reconstructions.py"]

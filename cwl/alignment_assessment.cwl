#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- run threaded alignment assessment
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)


inputs:
  assmb_map:
    label: Path to map from format_for_assembly.cwl
    type: File
    inputBinding:
      prefix: "--assmb_map"

  number_of_jobs:
    label: Number of alignment jobs to spawn
    type: int
    inputBinding:
      prefix: "--cpus"

  align_path:
    label: Path to output directory for all these alignments.
    type: Directory
    inputBinding:
      prefix: "--align_path"

  ivc_outfile:
    label: Path to the ids_v_cov.tsv outfile
    type: File
    inputBinding:
      prefix: "--ivc_outfile"

  best_only:
    label: Either "yes" or "no" for whether to report stats of only the best alignment or all alignments
    type: string
    inputBinding:
      prefix: "--best_only"

  assmb_type:
    label: Either "SPAdes" or "HGA". Determines how many assembled sequences are aligned to
    type: string
    inputBinding:
      prefix: "--assmb_type"

  priority:
    label: If given, the prefix of the sequence to solelys align to like XYZ.11203981.1 would require "XYZ" as input. Useful when trying to reconstruct a particular sequence
    type: string?
    inputBinding:
      prefix: "--priority"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?
    

outputs:
  ids_v_cov:
    type: File
    outputBinding:
      outputEval: $(inputs.ivc_outfile)


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/threaded_assess_alignment.py"]
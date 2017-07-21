#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- concatenate all the assembled sequence files
class: CommandLineTool


inputs:
  first_final:
    type: File
    inputBinding: {}

  second_final:
    type: File
    inputBinding: {}

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  all_seqs:
    type: stdout

stdout: 'assembled_seqs.fsa'


baseCommand: [cat]
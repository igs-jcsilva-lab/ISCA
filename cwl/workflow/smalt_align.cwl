#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: SMALT - Align to the index
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  smalt_dir:
    label: SMALT directory for the two stored files
    type: Directory

  prefix:
    label: SMALT prefix for the indices built in smalt_index.cwl
    type: string

  threads:
    label: number of threads to use for running SMALT
    type: int
    inputBinding:
      prefix: "-n"
      position: 1

  threshold:
    label: SMALT flag for deciding whether to keep suboptimal hits
    type: int
    inputBinding:
      prefix: "-d"
      position: 2
    default: -1

  smalt_indices:
    label: Path of the prefix for the indices built by smalt_index.cwl
    type: string
    inputBinding:
      position: 3
      valueFrom: $(inputs.smalt_dir.path + '/' + inputs.smalt_prefix)
    default: './smalt'

  reads1:
    label: Path to the first read pair file
    type: File
    inputBinding:
      position: 4

  reads2:
    label: Path to the second read pair file
    type: File
    inputBinding:
      position: 5

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  smalt_sam:
    type: stdout

stdout: 'smalt.sam'


baseCommand: ["/usr/local/packages/smalt/bin/smalt","map"]
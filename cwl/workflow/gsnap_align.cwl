#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: GSNAP - Align the paired reads to the index
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string
    inputBinding:
      prefix: "-d"
      position: 1

  gsnap_dir:
    label: Path to the output directory to write the GSNAP genome files to
    type: Directory
    inputBinding:
      prefix: "-D"
      position: 2

  threads:
    label: Number of threads to use for alignment
    type: int
    inputBinding:
      prefix: "-t"
      position: 3

  format:
    label: Format to output GSNAP alignment data
    type: string
    inputBinding:
      prefix: "-A"
      position: 4
    default: "sam"

  multihits:
    label: How many suboptimal hits to return
    type: int
    inputBinding:
      prefix: "-M"
      position: 5
    default: 6000

  reads1:
    label: Path to the first read pair file
    type: File
    inputBinding:
      prefix: "--gunzip"
      position: 6

  reads2:
    label: Path to the second read pair file
    type: File
    inputBinding:
      position: 7

  python3_lib:
    label: Path tp allow Python3 to be found in the ENV
    type: string?


outputs:
  gsnap_sam:
    type: stdout

stdout: 'gsnap.sam'


baseCommand: ["/usr/local/packages/gmap-gsnap/bin/gsnap"]
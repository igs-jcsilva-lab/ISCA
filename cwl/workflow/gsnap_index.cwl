#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: GSNAP - Build an index
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  gsnap_genome:
    inputBinding:
      prefix: "-d"
      position: 1
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string

  gsnap_dir:
    inputBinding:
      prefix: "-D"
      position: 2
    label: Path to the output directory to write the GSNAP genome files to
    type: string

  sequences:
    label: Path to the sequence file built from extract_sequences.py
    type: string
    inputBinding:
      position: 3


outputs:
  gsnap_out:
    type: Directory
    outputBinding:
      glob: $(inputs.gsnap_dir + '/' + inputs.gsnap_genome)

  gsnap_base_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.gsnap_dir)


baseCommand: ["/usr/local/packages/gmap-gsnap/bin/gmap_build"]
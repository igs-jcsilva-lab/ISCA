#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Establish workspace, first step of the pipeline
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  workspace_location:
    inputBinding:
      prefix: "-workspace_location"
    label: Path to build the workspace at, will write directories and pull new Python scripts here
    type: string
    default: '.'


outputs:
  gsnap_idx:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/gsnap_idx')

  smalt_idx:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/smalt_idx')

  sam:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/sam')

  reads:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/reads')

  spades_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/spades_assemblies')

  hga_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/hga_assemblies')

  alignments:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/alignments')

  second_reads:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/second_reads')

  second_spades_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/second_spades_assemblies')

  second_hga_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/second_hga_assemblies')

  second_alignments:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/second_alignments')

  HGA:
    type: File
    outputBinding:
      glob: $(inputs.workspace_location + '/HGA.py')

  scaffold_builder:
    type: File
    outputBinding:
      glob: $(inputs.workspace_location + '/scaffold_builder.py')


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/build_workspace.py"]
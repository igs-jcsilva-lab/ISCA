#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- Establish workspace, first step of the pipeline
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)


inputs:
  workspace_location:
    label: Path to build the workspace at, will write directories and pull new Python scripts here
    type: string
    inputBinding:
      prefix: "--workspace_location"
    default: "."
    
  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  gsnap_idx:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/gsnap_idx')

  smalt_idx:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/smalt_idx')

  sam_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/sam')

  first_reads:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/first_reads')

  first_spades_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/first_spades_assemblies')

  first_hga_assemblies:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/first_hga_assemblies')

  first_alignments:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/first_alignments')

  first_end_results:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/first_end_results')

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

  second_end_results:
    type: Directory
    outputBinding:
      glob: $(inputs.workspace_location + '/second_end_results')

  HGA:
    type: File
    outputBinding:
      glob: $(inputs.workspace_location + '/HGA.py')

  scaffold_builder:
    type: File
    outputBinding:
      glob: $(inputs.workspace_location + '/scaffold_builder.py')


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/build_workspace.py"]
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- run parallel assembly (SPAdes, HGA, ScaffoldBuilder)
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  python3_lib:
    label: Python3 library
    type: string

  reads_dir:
    label: Reads directory made by build_workspace.cwl
    type: Directory
    inputBinding:
      prefix: "--reads_dir"

  assmb_path:
    label: Path to the the directory to initialize directories for all the assembly output
    type: Directory
    inputBinding:
      prefix: "--assmb_path"

  assmb_map:
    label: Path to map from format_for_assembly.cwl
    type: File
    inputBinding:
      prefix: "--assmb_map"

  assmb_step:
    label: Either "gene" or "exon" for which sequences to pull
    type: string
    inputBinding:
      prefix: "--assmb_step"

  spades_install:
    label: Location of the SPAdes installation
    type: Directory
    inputBinding:
      prefix: "--spades_install"

  number_of_jobs:
    label: Number of assembly jobs to spawn
    type: int
    inputBinding:
      prefix: "--number_of_jobs"

  threads_per_job:
    label: Number of threads to use for each assembly job
    type: int
    inputBinding:
      prefix: "--threads_per_job"

  memory_per_job:
    label: How much memory to limit for each individual assembly job
    type: int
    inputBinding:
      prefix: "--memory_per_job"


outputs:
  assembled_dir:
    type: Directory
    outputBinding:
      outputEval: $(inputs.assmb_path)


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/run_parallel_assembly.py"]
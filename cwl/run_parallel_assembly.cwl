#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- run parallel assembly (SPAdes, HGA, ScaffoldBuilder)
class: CommandLineTool
stdout: run_parallel_assembly.out

requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  reads_dir:
    label: Reads directory made by build_workspace.cwl
    type: Directory?
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
    label: Either "SPAdes","HGA", or "SB"
    type: string
    inputBinding:
      prefix: "--assmb_step"

  spades_install:
    label: Location of the SPAdes installation
    type: Directory?
    inputBinding:
      prefix: "--spades_install"

  HGA_exe:
    label: Location of the HGA installation
    type: File?
    inputBinding:
      prefix: "--HGA_install"

  SB_exe:
    label: Location of the HGA installation
    type: File?
    inputBinding:
      prefix: "--SB_install"

  python2_exe:
    label: Location of the Python2 installation
    type: File?
    inputBinding:
      prefix: "--python2_install"

  velvet_install:
    label: Location of the Velvet installation
    type: Directory?
    inputBinding:
      prefix: "--velvet_install"

  number_of_jobs:
    label: Number of assembly jobs to spawn
    type: int
    inputBinding:
      prefix: "--number_of_jobs"

  partitions:
    label: Number of partitions to use in HGA
    type: int?
    inputBinding:
      prefix: "--partitions"

  threads_per_job:
    label: Number of threads to use for each assembly job
    type: int?
    inputBinding:
      prefix: "--threads_per_job"

  memory_per_job:
    label: How much memory to limit for each individual assembly job
    type: int?
    inputBinding:
      prefix: "--memory_per_job"

  ea_map:
    label: Path to the file built from extract_alleles.py
    type: File?
    inputBinding:
      prefix: "--ea_map"

  original_fsa:
    label: Path to where the unbuffered FASTA file generated from extract_sequences is
    type: File?
    inputBinding:
      prefix: "--original_fsa"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?

  assmb_stdout:
    label: Prior step stdout used as workaround for outputEval.
    type: File?

outputs:
  stdout:
    type: stdout


baseCommand: ["PYTHON3_EXE","TARGETED_ASSEMBLY_BIN/run_parallel_assembly.py"]
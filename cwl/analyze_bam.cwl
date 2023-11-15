#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- analyze the SAM/BAM file from aligner
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  sam:
    label: Path to a SAM file derived from aligner
    type: File
    inputBinding:
      prefix: "--sam"

  ea_map:
    label: Path to the output from extract_alleles.py
    type: File
    inputBinding:
      prefix: "--ea_map"

  threshold:
    label: Percent cut-off for how many matches a read must hit (uses whichever is longer the read or reference as denominator)
    type: int
    inputBinding:
      prefix: "--threshold"

  prefix:
    label: Name of the prefix to yield the two maps (one read-based and one reference-based)
    type: string
    inputBinding:
      prefix: "--prefix"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?

  samtools_install:
    label: Location of samtools
    type: Directory
    inputBinding:
      prefix: "--samtools_install"

outputs:
  read_map:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*read*")

  ref_map:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*ref*")


baseCommand: ["PYTHON3_EXE","TARGETED_ASSEMBLY_BIN/analyze_bam.py"]
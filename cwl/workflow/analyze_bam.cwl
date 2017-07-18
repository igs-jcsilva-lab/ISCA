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
  bam:
    label: Path to a BAM file derived from aligner
    type: File?
    inputBinding:
      prefix: "--bam"

  sam:
    label: Path to a SAM file derived from aligner
    type: File?
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


outputs:
  read_map:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*read*")

  ref_map:
    type: File
    outputBinding:
      glob: $(inputs.prefix + "*ref*")


baseCommand: ["/usr/local/packages/python-3.5.2/bin/python","/local/scratch/matsu_cwl_tests/analyze_bam.py"]
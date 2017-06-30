#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- analyze the SAM/BAM file from aligner
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


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


outputs:
  analyze_bam_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.prefix + "*")


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/analyze_bam.py"]
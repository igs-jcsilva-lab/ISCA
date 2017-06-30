#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- assign reads to individual locus directories
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
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

  ab_ref_map:
    label: Path to *_ref_map.tsv output from analyze_bam.cwl
    type: File
    inputBinding:
      prefix: "--ref_map"

  outfile:
    label: Name of the map to create which maps a reference locus to an int ID
    type: string
    inputBinding:
      prefix: "--outfile"


outputs:
  renamed_reads_dir:
    type: Directory
    outputBinding:
      outputEval: $(inputs.reads_dir)

  renamed_assmb_dir:
    type: Directory
    outputBinding:
      outputEval: $(inputs.assmb_path)

  assmb_map:
    type: File
    outputBinding:
      glob: $(inputs.outfile)


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/format_for_assembly.py"]
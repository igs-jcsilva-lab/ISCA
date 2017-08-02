#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- get the final assembled sequences
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)


inputs:
  ivc:
    label: Path to ids_v_cov map from assembly_verdict.cwl
    type: File
    inputBinding:
      prefix: "--ivc"

  threshold:
    label: Cutoff for pulling a sequence or not (% ID to reference)
    type: int
    inputBinding:
      prefix: "--threshold"

  align_path:
    label: Path to output directory for all these alignments.
    type: Directory
    inputBinding:
      prefix: "--align_path"

  groupby:
    label: Get sequences by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds"
    type: string
    inputBinding:
      prefix: "--groupby"

  ea_map:
    label: Path to the output from extract_alleles.py
    type: File
    inputBinding:
      prefix: "--ea_map"

  outfile:
    label: Name of output file
    type: string
    inputBinding:
      prefix: "--outfile"

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string?


outputs:
  final_sequences:
    type: File
    outputBinding:
      glob: $(inputs.outfile)


baseCommand: ["PYTHON3_EXE","TARGETED_ASSEMBLY_BIN/get_final_sequences.py"]
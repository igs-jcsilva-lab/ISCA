#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Phase three of the workflow consists of alignment.cwl, alignment_assessment.cwl, assembly_verdict.cwl, and get_final_sequences.cwl
class: Workflow


requirements:
  - class: InlineJavascriptRequirement


inputs:
  ea_map:
    label: Path to the output from extract_alleles.py
    type: File

  assmb_map:
    label: Path to map from format_for_assembly.cwl
    type: File

  ref_genome:
    label: Path to the reference genome file used to build aligner index
    type: File

  emboss_tool:
    label: Path to install directory of EMBOSS needle/water executable (e.g. /path/to/packages/emboss/bin/[needle|water])
    type: File

  original_fsa:
    label: Path to where the initial FASTA file generated from the pipeline is
    type: File

  assmb_path:
    label: Path to the the directory to initialize directories for all the assembly output
    type: Directory

  align_path:
    label: Path to output directory for all these alignments.
    type: Directory

  out_dir:
    label: Path to where the unaligned/unassembled FASTA entries and the new alignments map should go
    type: Directory

  assmb_type:
    label: Either "SPAdes" or "HGA". Determines how many assembled sequences are aligned to
    type: string

  ivc_outfile:
    label: Name of output file
    type: string

  sequences_outfile:
    label: Name of output file
    type: string

  best_only:
    label: Either "yes" or "no" for whether to report stats of only the best alignment or all alignments
    type: string

  groupby:
    label: Get sequences by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds"
    type: string

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string

  number_of_jobs:
    label: Number of alignment jobs to spawn
    type: int

  threshold:
    label: Minimum threshold of %ID that needs to be met to pass final assembly
    type: int

  min_align_len:
    label: Minimum alignment length to perform alignment with (RATIO)
    type: double


outputs:
  alignments:
    type: Directory
    outputSource: alignment/aligned_dir

  ids_v_cov:
    type: File
    outputSource: alignment_assessment/ids_v_cov

  final_results:
    type: Directory
    outputSource: assembly_verdict/end_results

  final_sequences:
    type: File
    outputSource: get_final_sequences/final_sequences


steps:
  alignment:
    run: alignment.cwl
    in:
      ea_map: ea_map
      assmb_map: assmb_map
      original_fsa: original_fsa
      assmb_path: assmb_path
      emboss_tool: emboss_tool
      align_path: align_path
      number_of_jobs: number_of_jobs
      min_align_len: min_align_len
      assmb_type: assmb_type
      python3_lib: python3_lib
    out: [aligned_dir]

  alignment_assessment:
    run: alignment_assessment.cwl
    in:
      assmb_map: assmb_map
      align_path: alignment/aligned_dir
      number_of_jobs: number_of_jobs
      assmb_type: assmb_type
      outfile: ivc_outfile
      best_only: best_only
      python3_lib: python3_lib
    out: [ids_v_cov]

  assembly_verdict:
    run: assembly_verdict.cwl
    in:
      ivc: alignment_assessment/ids_v_cov
      threshold: threshold
      original_fsa: original_fsa
      original_assmb_map: assmb_map
      out_dir: out_dir
      python3_lib: python3_lib
    out: [end_results]
    
  get_final_sequences:
    run: get_final_sequences.cwl
    in:
      ivc: alignment_assessment/ids_v_cov
      align_path: alignment/aligned_dir
      groupby: groupby
      threshold: threshold
      outfile: sequences_outfile
      ea_map: ea_map
      python3_lib: python3_lib
    out: [final_sequences]
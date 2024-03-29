#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Phase one of the workflow consists of build_workspace.cwl, extract_alleles.cwl, and extract_sequences.cwl
class: Workflow


requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  ea_input:
    label: Path to a TSV list for references and isolates
    type: File

  subset_list:
    label: Path to a file to subset sequences by
    type: File

  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int

  prefix:
    label: Prefix of the output FASTA file to generate in current or existing directory
    type: string

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string


outputs:
  gsnap_idx:
    type: Directory
    outputSource: build_workspace/gsnap_idx

  smalt_idx:
    type: Directory
    outputSource: build_workspace/smalt_idx

  first_reads:
    type: Directory
    outputSource: build_workspace/first_reads

  first_spades_assemblies:
    type: Directory
    outputSource: build_workspace/first_spades_assemblies

  first_hga_assemblies:
    type: Directory
    outputSource: build_workspace/first_hga_assemblies

  first_alignments:
    type: Directory
    outputSource: build_workspace/first_alignments

  second_reads:
    type: Directory
    outputSource: build_workspace/second_reads

  second_spades_assemblies:
    type: Directory
    outputSource: build_workspace/second_spades_assemblies

  second_hga_assemblies:
    type: Directory
    outputSource: build_workspace/second_hga_assemblies

  second_alignments:
    type: Directory
    outputSource: build_workspace/second_alignments

  first_ivc:
    type: File
    outputSource: build_workspace/first_ivc

  second_ivc:
    type: File
    outputSource: build_workspace/second_ivc

  ea_map:
    type: File
    outputSource: extract_alleles/ea_map

  buffered_sequences:
    type: File
    outputSource: extract_sequences/buffered_sequences

  unbuffered_sequences:
    type: File
    outputSource: extract_sequences/unbuffered_sequences


steps:
  build_workspace:
    run: build_workspace.cwl
    in:
      python3_lib: python3_lib
    out: [
      gsnap_idx,
      smalt_idx,
      first_reads,
      first_spades_assemblies,
      first_hga_assemblies,
      first_alignments,
      second_reads,
      second_spades_assemblies,
      second_hga_assemblies,
      second_alignments,
      first_ivc,
      second_ivc
    ]

  extract_alleles:
    run: extract_alleles.cwl
    in:
      ea_input: ea_input
      gene_or_exon: gene_or_exon
      python3_lib: python3_lib
    out: [ea_map]

  extract_sequences:
    run: extract_sequences.cwl
    in:
      ea_input: ea_input
      ea_map: extract_alleles/ea_map
      subset_list: subset_list
      buffer: buffer
      prefix: prefix
      python3_lib: python3_lib
    out: [
      buffered_sequences,
      unbuffered_sequences
    ]
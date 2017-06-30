#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Phase one of the workflow consists of build_workspace.cwl, extract_alleles.cwl, and extract_sequences.cwl
class: Workflow


requirements:
  - class: InlineJavascriptRequirement


inputs:
  workspace_location:
    label: Path to build the workspace at, will write directories and pull new Python scripts here
    type: string

  ea_input:
    label: Path to a TSV list for references and isolates
    type: File

  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int

  outfile:
    label: Name of the output FASTA file to generate in current or existing directory
    type: string


outputs:
  gsnap_idx:
    type: Directory
    outputSource: build_workspace/gsnap_idx

  smalt_idx:
    type: Directory
    outputSource: build_workspace/smalt_idx

  sam_dir:
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
  
  HGA:
    type: File
    outputSource: build_workspace/HGA
  
  scaffold_builder:
    type: File
    outputSource: build_workspace/scaffold_builder

  ea_map:
    type: File
    outputSource: extract_alleles/ea_map

  sequences:
    type: File
    outputSource: extract_sequences/sequences


steps:
  build_workspace:
    run: build_workspace.cwl
    in:
      workspace_location: workspace_location
    out: [
      gsnap_idx,
      smalt_idx,
      sam_dir,
      first_reads,
      first_spades_assemblies,
      first_hga_assemblies,
      first_alignments,
      second_reads,
      second_spades_assemblies,
      second_hga_assemblies,
      second_alignments,
      HGA,
      scaffold_builder
    ]

  extract_alleles:
    run: extract_alleles.cwl
    in:
      ea_input: ea_input
      gene_or_exon: gene_or_exon
    out: [ea_map]

  extract_sequences:
    run: extract_sequences.cwl
    in:
      ea_input: ea_input
      ea_map: extract_alleles/ea_map
      buffer: buffer
      outfile: outfile
    out: [sequences]
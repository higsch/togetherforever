#!/usr/bin/env nextflow

// set outdir
params.outdir = "results"

// read in gtf files from StringTie
// params.gtf has to be in parantheses
Channel
  .fromPath(params.gtfs)
  .map { it -> [it.baseName.split("\\.")[0], file(it)] }
  .set { gtfs }

// read genome fasta
genome_fasta = file(params.genome)

// set threeFrameTranslator.py
translator = file("threeFrameTranslator.py")


// fetch the nucleotide sequences from gtf file based on genome fasta
process getNucleotideSequences {

  input:
  set val(sample), file(gtf) from gtfs
  file genome_fasta

  output:
  set val("${sample}"), file("${sample}.fasta") into nucleotide_fastas

  script:
  """
  gffread -F -w ${sample}.fasta -g $genome_fasta ${gtf}
  """

}


// translate the nucleotide transcripts to amino acids (three frames)
process threeFrameTranslation {

  publishDir params.outdir, mode: "copy"
  
  input:
  set val(sample), file(nucleotide_fasta) from nucleotide_fastas
  file translator

  output:
  set val("${sample}"), file("${sample}.prot.fasta") into aa_fastas

  script:
  """
  python $translator -i $nucleotide_fasta -o ${sample}.prot.fasta
  """

}

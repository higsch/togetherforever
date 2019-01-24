#!/usr/bin/env nextflow

// set default outdir
params.outdir = "results"

// read in gtf files from StringTie and define sample names
// params.gtf has to be in parantheses
Channel
  .fromPath(params.gtfs)
  .map { it -> [it.baseName.split("\\.")[0], file(it)] }
  .set { gtfs }

// read genome fasta
genome_fasta = file(params.genome)

// canonical proteins
canonical_proteins = file(params.canonical)

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

aa_fastas
  .map { it -> it[1] }
  .collect()
  .set { aa_fastas_combined }

// merge all samples and remove duplicate IDs
process mergeSamplesFasta {

  input:
  file fastas from aa_fastas_combined

  output:
  file 'combined_unique.fasta' into fasta_combined_unique

  script:
  """
  for fasta in $fastas; do
    cat \${fasta} >> combined.fasta
  done
  awk 'NR%2 && !a[\$0]++ { print; getline l ; print l }' combined.fasta > combined_unique.fasta
  """

}

// add canonical proteins
// Todo: check awk command
process addCanonicalProteins {

  input:
  file 'combined_unique.fasta' from fasta_combined_unique
  file canonical_proteins

  output:
  file 'combined_unique_canonical.fasta' into fasta_combined_unique_canonical

  script:
  """
  cat $canonical_proteins combined_unique.fasta > tmp.fasta
  awk '((NR+1)%2) && !a[\$0]++ { print; getline l ; print l }' tmp.fasta > combined_unique_canonical.fasta
  """

}

// digest proteins
process digest {

  input:
  file 'combined_unique_canonical.fasta' from fasta_combined_unique_canonical

  output:
  file 'peptides.fasta' into peptides

  script:
  """
  /Applications/OpenMS-2.3.0/bin/Digestor -in combined_unique_canonical.fasta \
                                          -out peptides.fasta \
                                          -out_type fasta \
                                          -missed_cleavages 2 \
                                          -enzyme Trypsin/P
  """

}

// assign pI values by piDeepNet

// split database based on HiRIEF settings

// run MSGF+

// percolator

// protein inference

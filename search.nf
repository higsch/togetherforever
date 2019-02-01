#!/usr/bin/env nextflow

// default values
params.outdir = 'results'
params.missed_cleavages = 2
params.enzyme = 'Trypsin/P'

params.intercept = 3.5
params.width = 0.07
params.tolerance = 0.11
params.amount = 72

// read in gtf files from StringTie and define sample names
// params.gtf has to be in parantheses
Channel
  .fromPath(params.gtfs)
  .map { it -> [it.baseName.split('\\.')[0], file(it)] }
  .set { gtfs }

// read in mzML definition file
// layout: filepath, set, fraction
// params.mzmls has to be in parantheses
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize(' |\t') }
  .map { it -> [it[1], it[2], file(it[0])] } // set; fraction; file
  .tap { mzmls }
  .collect { it[1] }
  .set { fractions }

// read in normal PSMs from presearch
normPsms = file(params.normpsms)

// read genome fasta
genome_fasta = file(params.genome)

// canonical proteins
canonical_proteins_fasta = file(params.canonical)

// set executables
translator = file('threeFrameTranslator.py')
codonsplitter = file('codonsplitter.py')
piDeepNet = file('piDeepNet/getpiScores.R')
dbsplitter = file('dbsplitter.py')


// fetch the nucleotide sequences from gtf file based on genome fasta
process GetNucleotideSequences {

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
process ThreeFrameTranslation {
  
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
process MergeSampleFastas {

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


// split sequences with stop codon in separate sequences
process SplitStopCodons {

  input:
  file stop_fasta from fasta_combined_unique

  output:
  file 'no_stop.fasta' into fasta_nostop

  script:
  """
  python $codonsplitter -i $stop_fasta -o no_stop.fasta -c \"*\"
  """

}


// digest proteins
process DigestTranscriptome {

  input:
  file proteins_fasta from fasta_nostop

  output:
  file 'peptides.fasta' into peptides

  script:
  """
  Digestor -in $proteins_fasta \
           -out peptides.fasta \
           -out_type fasta \
           -missed_cleavages $params.missed_cleavages \
           -enzyme $params.enzyme
  """

}


// assign pI values by piDeepNet
// includes starting a h2o server before
process PIPredictionOnTranscriptome {

  input:
  file peptides from peptides

  output:
  file 'peptides_pI.fasta' into peptides_pI

  script:
  """
  java -jar /Users/matthias.stahl/ki/togetherforever/piDeepNet/piDeep/h2o.3.14.0.3/h2o/java/h2o.jar &
  Rscript $piDeepNet $peptides peptides_pI.fasta
  """

}


// split database based on isoelectric focussing
process SplitTranscriptomePeptidesToPIDBs {

  input:
  file peptides_pI from peptides_pI
  val fractions from fractions
  file normPsms

  output:
  file 'db_*' into pI_fastas

  script:
  """
  python $dbsplitter --pi-peptides $peptides_pI \
                     --normpsms $normPsms \
                     --intercept $params.intercept \
                     --width $params.width \
                     --tolerance $params.tolerance \
                     --amount $params.amount \
                     --fractions ${fractions.join(',')} \
                     --out db_*.fasta
  """

}


pI_fastas
  .flatten()
  .map { it -> [it.baseName.split("_")[1], file(it)] }
  .set { pI_tdbs }


// digest canonical proteins
process DigestKnownProteome {

  input:
  file canonical_proteins_fasta

  output:
  file 'canonical_peptides.fasta' into canonical_peptides

  script:
  """
  Digestor -in $canonical_proteins_fasta \
           -out canonical_peptides.fasta \
           -out_type fasta \
           -missed_cleavages $params.missed_cleavages \
           -enzyme $params.enzyme
  """

}


// add canonical peptides to each pI DB
process MergeTranscriptomeCanonicalsAndAddDecoys {

  publishDir 'results', mode: "copy"

  input:
  set val(fraction), file(db) from pI_tdbs
  file(canonical_peptides) from canonical_peptides

  output:
  set val("${fraction}"), file("tddb_${fraction}.fasta") into combined_tdbs

  script:
  """
  DecoyDatabase -in $db $canonical_peptides \
                -out tddb_${fraction}.fasta
  """

}

// run MSGF+

// percolator

// protein inference

#!/usr/bin/env nextflow

// default values
params.outdir = "results"
params.missed_cleavages = 2
params.enzyme = "Trypsin/P"

// read in gtf files from StringTie and define sample names
// params.gtf has to be in parantheses
Channel
  .fromPath(params.gtfs)
  .map { it -> [it.baseName.split("\\.")[0], file(it)] }
  .set { gtfs }

// read in mzML definition file
// layout: filepath, set, fraction
// params.mzmls has to be in parantheses
Channel
  .from(file(params.mzmldef).readLines()) // include $, if it doesn't work
  .map {it -> it.tokenize("\t") }
  .map { it -> [it[1], it[2], file(it[0])] } // set; fraction; file
  .set { mzmls }

// read in normal PSMs from presearch
Channel
  .from(file(params.normalpsms))
  .set { normalPsms}

// read genome fasta
genome_fasta = file(params.genome)

// canonical proteins
canonical_proteins = file(params.canonical)

// set executables
translator = file("threeFrameTranslator.py")
codonsplitter = file("codonsplitter.py")
piDeepNet = file("piDeepNet/getpiScores.R")
dbsplitter = file("dbsplitter.py")


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

// add canonical proteins and delete duplicates
// adding the canonical to the top ensures that
// transcript-based duplicates are omitted
// test this later, without duplicate removal
process addCanonicalProteins {

  input:
  file 'combined_unique.fasta' from fasta_combined_unique
  file canonical_proteins

  output:
  file 'combined_unique_canonical.fasta' into fasta_combined_unique_canonical

  script:
  """
  # combine fastas
  cat $canonical_proteins combined_unique.fasta > tmp.fasta

  # make one line per header and sequence
  seqtk seq -A  < tmp.fasta > tmp_layouted.fasta

  # Check, if every second line is a header
  number_headers=`awk '(NR+1)%2==0' tmp_layouted.fasta | grep -e '^>' -c`
  number_lines=`wc -l < tmp_layouted.fasta`
  number_lines_half=\$(( number_lines / 2 ))

  # if true, remove sequence duplicates
  if [ \"\$number_headers\" = \"\$number_lines_half\" ]; then
    awk -v RS='>' -v ORS="" 'NR == 1 { next } seen[\$NF] { delete seen[\$NF]; next } { print ">" \$0; seen[\$NF]=1 }' tmp_layouted.fasta > combined_unique_canonical.fasta
    # awk 'BEGIN{RS=">";ORS="";} NF>0 && !a[\$NF]++ { print ">"\$0; }' tmp_layouted.fasta > combined_unique_canonical.fasta
  else
    echo \"Error! The input fastas are assumed to have one line per header and sequence.\"
  fi
  """

}


// split sequences with stop codon in separate sequences
process splitStopCodons {

  input:
  file 'combined_unique_canonical.fasta' from fasta_combined_unique_canonical

  output:
  file 'nostop.fasta' into fasta_nostop

  script:
  """
  python $codonsplitter -i combined_unique_canonical.fasta -o nostop.fasta -c \"*\"
  """

}


// digest proteins
process digestProteins {

  input:
  file 'nostop.fasta' from fasta_nostop

  output:
  file 'peptides.fasta' into peptides

  script:
  """
  Digestor -in nostop.fasta \
           -out peptides.fasta \
           -out_type fasta \
           -missed_cleavages $params.missed_cleavages \
           -enzyme $params.enzyme
  """

}


// assign pI values by piDeepNet
// start h2o server before
process piDeepNet {

  publishDir params.outdir, mode: "copy"

  input:
  file 'peptides.fasta' from peptides

  output:
  file 'peptides_pI.fasta' into peptides_pI

  script:
  """
  Rscript $piDeepNet peptides.fasta peptides_pI.fasta
  """

}


// split database based on HiRIEF settings
process splitPeptidesToPIFastas {

  input:
  file 'peptides_pI.fasta' from peptides_pI
  file normalPsms from normalPsms

  output:
  file 'db_*' into pI_fastas

  script:
  """
  python $dbsplitter 
  """

}

// run MSGF+

// percolator

// protein inference

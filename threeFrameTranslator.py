#!/usr/bin/env python
# from https://github.com/husensofteng/ProteoGenomics
import sys
import os
import glob
import argparse
from Bio import SeqIO

def parseCommandlineArgs():
  parser = argparse.ArgumentParser(description = '3-frame translation for nucleotide sequences.')
  parser.add_argument('-i', '--input-file', dest = "infile", default = "", required = True, help = 'input fasta')
  parser.add_argument('-o', '--output-file', dest = "outfile", default = "", required = True, help = 'output fasta')
  return parser.parse_args(sys.argv[1:])

def translate(infile, outfile):
  with open(infile, 'r') as handle, open(outfile, 'w') as output:
    for record in SeqIO.parse(handle, 'fasta'):
      seq = record.seq
      ORF1 = seq.translate()
      ORF2 = seq[1::].translate()
      ORF3 = seq[2::].translate()
      if record.id == "":
        print(record.description, record.seq)
        continue
      output.write("%s\n%s\n" % ('>'+record.id+"_ORF1 "+record.description, ORF1))
      output.write("%s\n%s\n" % ('>'+record.id+"_ORF2 "+record.description, ORF2))
      output.write("%s\n%s\n" % ('>'+record.id+"_ORF3 "+record.description, ORF3))

if __name__ == '__main__':
  args = parseCommandlineArgs()
  translate(args.infile, args.outfile)

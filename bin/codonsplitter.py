#!/usr/bin/env python
# from https://github.com/husensofteng/ProteoGenomics
import sys
import os
import glob
import argparse
from Bio import SeqIO


def parseCommandlineArgs():
  parser = argparse.ArgumentParser(description = '3-frame translation for nucleotide sequences.')
  parser.add_argument('-i', '--input-file', dest = 'infile', default = '', required = True, help = 'input fasta')
  parser.add_argument('-o', '--output-file', dest = 'outfile', default = '', required = True, help = 'output fasta')
  parser.add_argument('-c', '--codon', dest = 'codon', default = '', required = False, help = 'codon')
  return parser.parse_args(sys.argv[1:])


def split(infile, outfile, codon):
  with open(infile, 'r') as handle, open(outfile, 'w') as output:
    for record in SeqIO.parse(handle, 'fasta'):
      seq = record.seq
      split_seq = seq.split(codon)
      splitter_counter = 0
      for splitter in split_seq:
        if len(splitter) > 0:
          output.write('%s\n%s\n' % ('>' + record.id + '_split' + str(splitter_counter) + ' ' + record.description + " split=" + str(splitter_counter), splitter))
          splitter_counter += 1


if __name__ == '__main__':
  args = parseCommandlineArgs()
  split(args.infile, args.outfile, args.codon)

#!/usr/bin/env python
# from https://github.com/husensofteng/ProteoGenomics
import sys
import os
import glob
import argparse
from Bio import SeqIO

def get_abspath(p):
  if not os.path.exists(os.path.abspath(os.path.expanduser(p))):
    if len(glob.glob(p))>0:
      return p
    raise argparse.ArgumentTypeError(p)
  if os.path.isdir(os.path.abspath(os.path.expanduser(p))):
    return os.path.abspath(os.path.expanduser(p))+'/'
  return os.path.abspath(os.path.expanduser(p))

def parseCommandlineArgs():
  parser = argparse.ArgumentParser(description = '3-frame translation for nucleotide sequences.')
  parser.add_argument('-i', '--input_file', type = get_abspath, default = "", required = True, help = 'input fasta')
  parser.add_argument('-o', '--output_file', type = get_abspath, default = "", required = True, help = 'output fasta')
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
  translate(args.input_file, args.output_file)

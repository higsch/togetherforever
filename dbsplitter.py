#!/usr/bin/env python
import sys
import argparse
from numpy import median
from Bio import SeqIO
import re


def parseCommandlineArgs():
  parser = argparse.ArgumentParser(description = 'Database splitter based on pI.')
  parser.add_argument('--pi-peptides', dest = 'piPeptides', default = '', required = True, help = 'peptides with pI')
  parser.add_argument('--normpsms', dest = 'normPsms', default = '', required = True, help = 'psms from presearch')
  parser.add_argument('--intercept', dest = 'intercept', default = '', required = True, help = 'strip intercept')
  parser.add_argument('--width', dest = 'width', default = '', required = True, help = 'strip width')
  parser.add_argument('--tolerance', dest = 'tolerance', default = '', required = True, help = 'strip tolerance')
  parser.add_argument('--amount', dest = 'amount', default = '', required = True, help = 'amount')
  parser.add_argument('--fractions', dest = 'fractions', default = '', required = True, help = 'real fractions')
  parser.add_argument('--out', dest = 'outfile', default = '', required = True, help = 'out file')
  return parser.parse_args(sys.argv[1:])


# Todo: Use pandas to calculate median
def getPIShift(normpsms):
  deltaPIs = []
  with open(normpsms) as fp:
    next(fp)  # skip header
    for line in fp:
      line = line.strip('\n').split('\t')
      try:
        deltaPI = float(line[43])
      except ValueError:
        continue
      if deltaPI < 0.2:
        deltaPIs.append(deltaPI)
    shift = median(deltaPIs)
    print('pI shift (median of delta pIs): {}'.format(shift))
    return shift


def calculatePIRange(lowerFracNr, upperFracNr, intercept, width, tolerance, pIShift):
  pICenter = width * ((lowerFracNr + upperFracNr) / 2) + intercept - pIShift
  pILeft = (width * lowerFracNr) + intercept - (width / 2) - tolerance - pIShift
  pIRight = (width * upperFracNr) + intercept + (width / 2) + tolerance - pIShift
  return (pILeft, pICenter, pIRight)


def getPiFracDict(fractions, amount, intercept, width, tolerance, pIShift):
  pIFracDict = {}
  fractionsArr = fractions.split(',')[::-1]
  lastFrac = amount
  for frac in fractionsArr:
    pIRange = calculatePIRange(frac, lastFrac - 1, intercept, width, tolerance, pIShift)
    pIFracDict[frac] = {
      'leftPI': pIRange[0],
      'centerPI': pIRange[1],
      'rightPI': pIRange[2],
      'peptides': []
    }
    lastFrac = frac
  return pIFracDict


def getPIFromFastaDescription(header):
  pI = -1
  pIString = re.match(r'pI=(?.+)[\n\s]', header).strip()
  try:
    pI = float(pIString)
  except ValueError:
    print('pI could not be read from Fasta header!')
  return pI


def searchPIFrac(pIFracDict, pI):
  counter = 0
  for frac, db in pIFracDict.items:
    if (pI >= db['leftPI'] and pI <= db['rightPI']):
      if (counter >= 3): break
      counter += 1
      yield frac


def assignPeptidesToFracDict(pIFracDict, pIPeptides):
  with open(pIPeptides, 'r') as peps:
    for record in SeqIO.parse(peps, 'fasta'):
      pI = getPIFromFastaDescription(record.description)
      if (pI != -1):
        for frac in searchPIFrac(pIFracDict, pI):
          pIFracDict[frac]['peptides'].append(record)
  return pIFracDict


def writeFasta(pIFracDict, outfile):
  pass


if __name__ == '__main__':
  args = parseCommandlineArgs()
  pIShift = getPIShift(args.normpsms)
  pIFracDict = getPiFracDict(args.fractions, args.amount, args.intercept, args.width, args.tolerance, pIShift)
  pIFracDict = assignPeptidesToFracDict(pIFracDict, args.pIPeptides)
  writeFasta(pIFracDict, args.outfile)

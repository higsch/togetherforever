#!/usr/bin/env python
import sys
import argparse
from numpy import median
from Bio import SeqIO
import re


def parseCommandlineArgs():
  parser = argparse.ArgumentParser(description = 'Database splitter based on pI.')
  parser.add_argument('--pi-peptides', dest = 'pIPeptides', default = '', required = True, help = 'peptides with pI')
  parser.add_argument('--normpsms', dest = 'normPsms', default = '', required = True, help = 'psms from presearch')
  parser.add_argument('--intercept', dest = 'intercept', type = float, default = '', required = True, help = 'strip intercept')
  parser.add_argument('--width', dest = 'width', type = float, default = '', required = True, help = 'strip width')
  parser.add_argument('--tolerance', dest = 'tolerance', type = float, default = '', required = True, help = 'strip tolerance')
  parser.add_argument('--amount', dest = 'amount', type = int, default = '', required = True, help = 'amount')
  parser.add_argument('--fractions', dest = 'fractions', default = '', required = True, help = 'real fractions')
  parser.add_argument('--out', dest = 'outFile', default = '', required = True, help = 'out file')
  return parser.parse_args(sys.argv[1:])


# Todo: Use pandas to calculate median
def getPIShift(normpsms):
  deltaPIs = []
  with open(normpsms) as fp:
    next(fp)  # skip header
    for line in fp:
      line = line.strip('\n').split('\t')
      try:
        deltaPI = float(line[42])
      except ValueError:
        continue
      if deltaPI < 0.2:
        deltaPIs.append(deltaPI)
    shift = median(deltaPIs)
    # print('pI shift (median of delta pIs): {}'.format(shift))
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
    frac = int(frac)
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
  pIString = re.search('(?<=pI=)\d+\\.\d+', header)
  try:
    pI = float(pIString.group(0))
  except ValueError:
    print('pI could not be read from Fasta header!')
  return pI


def searchPIFrac(pIFracDict, pI):
  fracs = []
  for frac, db in pIFracDict.items():
    if (pI < db['leftPI']):
      continue
    elif (pI <= db['rightPI']):
      fracs.append(frac)
    else:
      return fracs
    return fracs


def assignPeptidesToFracDict(pIFracDict, pIPeptides):
  with open(pIPeptides, 'r') as peps:
    for record in SeqIO.parse(peps, 'fasta'):
      pI = getPIFromFastaDescription(record.description)
      if (pI != -1):
        fracs = searchPIFrac(pIFracDict, pI)
        if (fracs):
          for frac in fracs:
            pIFracDict[frac]['peptides'].append(record)
  return pIFracDict


def writeFasta(pIFracDict, outFile):
  outFileComponents = outFile.split('*')
  for frac, db in pIFracDict.items():
    with open(outFileComponents[0] + str(frac) + outFileComponents[1], 'w') as o:
      for record in db['peptides']:
        o.write('%s\n%s\n' % ('>' + record.description, record.seq))


if __name__ == '__main__':
  args = parseCommandlineArgs()
  pIShift = getPIShift(args.normPsms)
  pIFracDict = getPiFracDict(args.fractions, args.amount, args.intercept, args.width, args.tolerance, pIShift)
  pIFracDict = assignPeptidesToFracDict(pIFracDict, args.pIPeptides)
  writeFasta(pIFracDict, args.outFile)

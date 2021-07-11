import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import logging
import argparse
import os
import csv
import re
from pprint import pprint

def getProviralFastaIDs(fafile, recordSeqs):
  ids = []
  for record in SeqIO.parse(fafile, format = "fasta"):
    ids.append(record.id)
    recordSeqs[record.id].append(record.seq)

  return ids


def extractCellBarcode(read):
  # accept only CB tag because it passes the allowlist set by 10X
  tags = dict(read.tags)
  
  if "CB" in tags:
      barcode = tags["CB"]
  else:
      barcode = "None"

  return barcode


def getLTRseq(seq, start, end):
  ltrSeq = seq[start - 1:end]
  return ltrSeq


def parseLTRMatches(LTRmatchesFN, proviralSeqs, endBuffer = 20):
  LTRdict = defaultdict(lambda: {"5p": [], "5pRevComp": [], "3p": [], "3pRevComp":[]})

  with open(LTRmatchesFN, "r") as fhandle:
    rd = csv.reader(fhandle, delimiter = "\t")

    for row in rd:
      # index 1 = subject ID (i.e. the original sample's viral fasta ID)
      # index 2 = percent match
      # index 6 = query start
      # index 7 = query end
      # index 8 = subject start
      # index 9 = subject end
      
      print(row)

      subjID = row[1]
      # qstart = row[6]
      # qend = row[7]
      sstart = int(row[8])
      send = int(row[9])
      slen = len(proviralSeqs[subjID][0])
      
      print(subjID, sstart, send)

      if sstart < endBuffer:
        seq = getLTRseq(proviralSeqs[subjID][0], 1, send)
        LTRdict[subjID]["5p"].append(seq)
        LTRdict[subjID]["5pRevComp"].append(seq.reverse_complement())
      elif slen - send < endBuffer:
        seq = getLTRseq(proviralSeqs[subjID][0], sstart, slen)
        LTRdict[subjID]["3p"].append(seq)
        LTRdict[subjID]["3pRevComp"].append(seq.reverse_complement())

  return LTRdict


def getSoftClip(read, clipMinLen):
  # cutoff same as epiVIA
  cigar = read.cigartuples
  clippedFrag = ""

  # clip is at 5' position
  if cigar[0][0] == 4 and cigar[0][1] >= clipMinLen:
    clipLen = cigar[0][1]
    clippedFrag = read.seq[0:clipLen]
  
  # clip at 3' position
  elif cigar[-1][0] == 4 and cigar[-1][0] >= clipMinLen:
    clipLen = cigar[-1][1]
    clippedFrag = read.seq[(cigar[-1][1]*-1): ]

  return clippedFrag


def isSoftClipProviral(read, proviralLTRSeqs, clipMinLen = 11):
  clippedFrag = getSoftClip(read, clipMinLen)
  
  if len(clippedFrag) == 0:
    return False

  hits = {"plus": [], "plusIds": [], "minus" : [], "minusIds": []}
  for key in proviralLTRSeqs:
    keyPair = proviralLTRSeqs[key]

    for ltrType in keyPair:
      orient = ""
      if ltrType == "5p" or ltrType == "3p":
        orient = "plus"
      else:
        orient = "minus"

      # find hits...
      ltrTypeSeqs = keyPair[ltrType]
      for s in ltrTypeSeqs:
        hits = [x.start() for x in re.finditer(clippedFrag, s)]
        if len(hits) > 0:
          hits[orient].append(hits)
          hits[orient + "Ids"].append(key)

    # check hits
    pprint(hits)
  
  return


def parseHostReadsWithPotentialChimera(reads, clipMinLen):
  for r in reads:
    isSoftClipProviral(r, clipMinLen)



def parseCellrangerBam(bamfile, proviralFastaIds, proviralReads, hostReadsWithPotentialChimera, unmappedPotentialChimera, top_n = -1):
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 20)

  readIndex = 0
  for read in bam:
    refnameIsProviral = read.reference_name in proviralFastaIds
    # supposed to take mate's ref name or if no mate, the next record in BAM file
    nextRefnameIsProviral = read.next_reference_name in proviralFastaIds
    
    cigarString = read.cigartuples
    # 4 is soft clip
    hasSoftClipAtEnd = cigarString[-1][0] == 4 or cigarString[0][0] == 4
    softClipInitThresh = 7
    softClipIsLongEnough = cigarString[-1][1] >= softClipInitThresh or cigarString[0][1] >= softClipInitThresh

    # ignore if optical/PCR duplicate OR without a mate
    if (read.flag & 1024) or (not read.flag & 1):
      readIndex += 1
      continue
    
    # if read is properly mapped in a pair AND not proviral aligned AND there is soft clipping involved
    # TODO add only if S number is > than settings
    elif (read.flag & 2) and (not refnameIsProviral) and (hasSoftClipAtEnd and softClipIsLongEnough):
      # move to chimera identification
      hostReadsWithPotentialChimera[read.qname].append(read)
    
    # if there is a mate AND both are proviral only 
    elif refnameIsProviral and nextRefnameIsProviral:
      # save into proviral
      proviralReads[read.qname].append(read)

    # read or mate must be mapped AND either read or its mate must be proviral
    elif (not read.flag & 14) and (refnameIsProviral or nextRefnameIsProviral):
      # move to chimera identification
      unmappedPotentialChimera[read.qname].append(read)
    
    readIndex += 1
    
    if readIndex % 1000000 == 0:
      logging.info("Parsed %s reads", str(readIndex))

    if top_n != -1 and readIndex > top_n:
      return
    
  return bam


def writeBam(fn, templateBam, reads):
  outputBam = pysam.AlignmentFile(fn, "wb", template = templateBam)

  for qname in reads:
    for read in reads[qname]:
      outputBam.write(read)


def main(args):
  # output filenames
  outputFNs = {
    "proviralReads": "proviralReads.bam",
    "hostWithPotentialChimera": "hostWithPotentialChimera.bam",
    "umappedWithPotentialChimera": "unmappedWithPotentialChimera.bam"
  }

  # set up initial dictionaries
  dualProviralAlignedReads = defaultdict(list)
  hostReadsWithPotentialChimera = defaultdict(list)
  unmappedPotentialChimera = defaultdict(list)

  # set up logger
  logging.basicConfig(filename='haystack.log', encoding='utf-8', level=logging.DEBUG)

  # recover all proviral "chromosome" names from partial fasta file used by Cellranger
  logging.info("Getting proviral records")
  proviralSeqs = defaultdict(lambda: [])
  proviralFastaIds = getProviralFastaIDs(args.viralFasta, proviralSeqs)

  # get possible LTR regions from fasta file
  # TODO build into program (as opposed to running it before)
  logging.info("Getting potential LTRs")
  potentialLTR = parseLTRMatches(args.LTRmatches, proviralSeqs)
  pprint(dict(potentialLTR))

  # parse BAM file
  logging.info("Parsing cellranger BAM (namesorted)")
  parseCellrangerBam(bamfile = args.bamfile,
    proviralFastaIds = proviralFastaIds,
    proviralReads = dualProviralAlignedReads,
    hostReadsWithPotentialChimera = hostReadsWithPotentialChimera,
    unmappedPotentialChimera = unmappedPotentialChimera,
    top_n = args.topNReads) #debugging

  # output BAM files
  logging.info("Writing out BAM files of parsed records")

  cellrangerBam = pysam.AlignmentFile(args.bamfile, "rb")
  writeBam(args.outputDir + "/" + outputFNs["proviralReads"], cellrangerBam, dualProviralAlignedReads)
  writeBam(args.outputDir + "/" + outputFNs["hostWithPotentialChimera"], cellrangerBam, hostReadsWithPotentialChimera)
  writeBam(args.outputDir + "/" +  outputFNs["umappedWithPotentialChimera"], cellrangerBam, unmappedPotentialChimera)
  cellrangerBam.close()

  # parse host reads with potential chimera
  parseHostReadsWithPotentialChimera(hostReadsWithPotentialChimera, clipMinLen == args.LTRClipLen)


if __name__ == '__main__':
  # set up command line arguments
  parser = argparse.ArgumentParser(
    description = "Identify HIV-associated reads amidst cellular scATAC reads")

  parser.add_argument("--bamfile",
    required = True,
    help = "Name sorted Cellranger BAM file")
  parser.add_argument("--outputDir",
    required = True,
    help = "Output bam files")
  parser.add_argument("--viralFasta",
    required = True,
    help = "Viral fasta file (can have multiple sequences in single file)")
  parser.add_argument("--topNReads",
    default = -1,
    type = int,
    help = "Limit to n number of records in BAM file. Default is all (-1)")
  parser.add_argument("--LTRmatches",
    required = True,
    help = "blastn table output format for LTR matches to HXB2 LTR")
  parser.add_argument("--LTRClipLen",
    default = 11,
    type = int,
    help = "Number of bp to extend into LTR from a chimeric fragment")
  args = parser.parse_args()

  if not os.path.exists(args.outputDir):
    os.makedirs(args.outputDir)

  main(args)

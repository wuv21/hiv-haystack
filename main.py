import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
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
      
      subjID = row[1]
      # qstart = row[6]
      # qend = row[7]
      sstart = int(row[8])
      send = int(row[9])
      slen = len(proviralSeqs[subjID][0])
      
      # must be at least 550bp long
      if abs(send - sstart) < 550:
        continue
      elif sstart < endBuffer:
        seq = getLTRseq(proviralSeqs[subjID][0], 1, send)
        LTRdict[subjID]["5p"].append(seq)
        LTRdict[subjID]["5pRevComp"].append(seq.reverse_complement())
      elif slen - send < endBuffer:
        seq = getLTRseq(proviralSeqs[subjID][0], sstart, slen)
        LTRdict[subjID]["3p"].append(seq)
        LTRdict[subjID]["3pRevComp"].append(seq.reverse_complement())

  return LTRdict


def getSoftClip(read, clipMinLen, softClipPad):
  # cutoff same as epiVIA
  cigar = read.cigartuples
  clippedFrag = Seq("")
  adjacentFrag = Seq("")

  clip5Present = False
  clip3Present = False

  # loop through cigar to make sure there's only 1 soft clip
  softClipCount = 0
  for c in cigar:
    if c[0] == 4:
      softClipCount += 1

  if softClipCount > 1:
    return None

  # clip is at 5' position
  if cigar[0][0] == 4 and cigar[0][1] >= clipMinLen:
    clipLen = cigar[0][1]
    clippedFrag = read.seq[0:clipLen]
    adjacentFrag = read.seq[clipLen:clipLen + softClipPad]
    clip5Present = True

  # clip at 3' position
  if cigar[-1][0] == 4 and cigar[-1][1] >= clipMinLen:
    clipLen = cigar[-1][1]
    clippedFrag = read.seq[clipLen * -1: ]
    adjacentFrag = read.seq[clipLen - softClipPad : clipLen * -1]
    clip3Present = True

  # clip can only be present at one end
  if clip5Present and clip3Present:
    return None
  else:
    clippedFragObj = {
      "clippedFrag": clippedFrag,
      "adjacentFrag": adjacentFrag,
      "clip5Present": clip5Present,
      "clip3Present": clip3Present}

    return clippedFragObj


def isSoftClipProviral(read, proviralLTRSeqs, clipMinLen = 11, softClipPad = 3):
  clippedFragObj = getSoftClip(read, clipMinLen, softClipPad)
  
  # skip if no clipped fragment long enough is found
  if clippedFragObj is None:
    return False

  # check if mate is in correct orientation
  # This actually should already be ok because the read pair already has to be properly mapped.
  # if clippedFragObj["clip5Present"] and read.next_reference_start < read.reference_start:
  #   return False
  # elif clippedFragObj["clip3Present"] and read.next_reference_start > read.reference_start:
  #   return False

  strClippedFrag = str(clippedFragObj["clippedFrag"])

  # skip if there are any characters other than ATGC 
  if bool(re.compile(r'[^ATGC]').search(strClippedFrag)):
    return False
  
  hits = {
    "plus": [],
    "plusIds": [],
    "minus" : [],
    "minusIds": [],
    "clip5P": clippedFragObj["clip5Present"],
    "clip3P": clippedFragObj["clip3Present"]}

  allowedLTRKeys = []
  # only allow specific keys based on orientation
  if clippedFragObj["clip5Present"]:
    allowedLTRKeys = ["3p", "5pRevComp"]
  elif clippedFragObj["clip3Present"]:
    allowedLTRKeys = ["5p", "3pRevComp"]

  # find hits...
  foundHit = False
  for key in proviralLTRSeqs:
    keyPair = proviralLTRSeqs[key]

    for ltrType in allowedLTRKeys:
      ltrTypeSeqs = keyPair[ltrType]
      if len(ltrTypeSeqs) == 0:
        continue

      # find orientation
      orient = "plus" if ltrType == "5p" or ltrType == "3p" else "minus"

      for s in ltrTypeSeqs:
        matches = [x.start() for x in re.finditer(strClippedFrag, str(s))]
        if len(matches) == 0:
          continue

        ltrLen = len(str(s))
        # check if match is within soft buffer zone
        # needs to pass min(matches) <= softClipPad or max(matches) + len(strClippedFrag) >= ltrLen - softClipPad:
        if (ltrType == "5p" or ltrType == "3pRevComp") and min(matches) > softClipPad:
          continue
        elif (ltrType == "3p" or ltrType == "5pRevComp") and max(matches) + len(strClippedFrag) < ltrLen - softClipPad:
          continue

        # check if the adjacent host clips could have also been aligned to the viral LTR,
        # thus explaining the lack of viral clip not being at either end of LTR
        ltrEnd = ""
        if (ltrType == "5p" or ltrType == "3pRevComp") and min(matches) != 0:
          ltrEnd = str(s)[0:min(matches)]
          hostAdjacentBp = str(read.seq)[-len(strClippedFrag) - len(ltrEnd): -len(strClippedFrag)]

        elif (ltrType == "3p" or ltrType == "5pRevComp") and max(matches) != ltrLen - softClipPad:
          adjacentBpNum = ltrLen - max(matches) - len(strClippedFrag)
          ltrEnd = str(s)[max(matches) + len(strClippedFrag): ltrLen]
          hostAdjacentBp = str(read.seq)[len(strClippedFrag): len(strClippedFrag) + adjacentBpNum]

        if ltrEnd != "" and ltrEnd != hostAdjacentBp:
          print("Viral clip not found at the end of LTR")
          print(s, matches, read.seq, read.query_name)
          continue

        # passes all checks!
        print("successful match")
        print(s, matches, read.seq, read.query_name)
        hits[orient].append(matches)
        hits[orient + "Ids"].append(key + "___" + ltrType)
        foundHit = True

  # can only be plus orientation OR minus orientation only
  if foundHit and len(hits["plus"]) != 0 and len(hits["minus"]) == 0:
    return hits
  elif foundHit and len(hits["minus"]) != 0 and len(hits["plus"]) == 0:
    return hits
  else:
    return False


def parseHostReadsWithPotentialChimera(readPairs, proviralLTRSeqs, clipMinLen):
  validHits = []
  validReads = []

  for key in readPairs:
    # only allow one read mate to have soft clip
    if len(readPairs[key]) != 1:
      continue 
    
    read = readPairs[key][0]
    potentialHits = isSoftClipProviral(read, proviralLTRSeqs, clipMinLen)
    
    if potentialHits and len(potentialHits['plus']) != 0:
      validHits.append(potentialHits)
      validReads.append(read)
    elif potentialHits and len(potentialHits['minus']) != 0:
      validHits.append(potentialHits)
      validReads.append(read)

  returnVal = {"validHits": validHits, "validReads": validReads}
  return returnVal


def parseCellrangerBam(bamfile, proviralFastaIds, proviralReads, hostReadsWithPotentialChimera, unmappedPotentialChimera, top_n = -1):
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 20)

  readIndex = 0
  for read in bam:
    # ignore if optical/PCR duplicate OR without a mate
    if (read.flag & 1024) or (not read.flag & 1):
      readIndex += 1
      continue
    
    refnameIsProviral = read.reference_name in proviralFastaIds
    # supposed to take mate's ref name or if no mate, the next record in BAM file
    nextRefnameIsProviral = read.next_reference_name in proviralFastaIds
    
    cigarString = read.cigartuples
    # 4 is soft clip
    hasSoftClipAtEnd = cigarString != None and (cigarString[-1][0] == 4 or cigarString[0][0] == 4)
    softClipInitThresh = 7
    softClipIsLongEnough = cigarString != None and (cigarString[-1][1] >= softClipInitThresh or cigarString[0][1] >= softClipInitThresh)
    
    # if read is properly mapped in a pair AND not proviral aligned AND there is soft clipping involved
    if (read.flag & 2) and (not refnameIsProviral) and (hasSoftClipAtEnd and softClipIsLongEnough):
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
      print("Parsed {} reads".format(str(readIndex)))

    if top_n != -1 and readIndex > top_n:
      return
    
  return bam


def writeBam(fn, templateBam, reads):
  outputBam = pysam.AlignmentFile(fn, "wb", template = templateBam)
  
  if isinstance(reads, list):
    for read in reads:
      outputBam.write(read)
  else:
    for qname in reads:
      for read in reads[qname]:
        outputBam.write(read)


def importProcessedBam(bamfile, returnDict = True):
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 20)

  if returnDict:
    val = defaultdict(list)
  else:
    val = []

  for read in bam:
    if returnDict:
      val[read.qname].append(read)
    else:
      val.append(read)

  return val


def main(args):
  # output filenames
  outputFNs = {
    "proviralReads": "proviralReads.bam",
    "hostWithPotentialChimera": "hostWithPotentialChimera.bam",
    "umappedWithPotentialChimera": "unmappedWithPotentialChimera.bam",
    "hostWithValidChimera": "hostWithValidChimera.bam"
  }
  
  # set up initial dictionaries
  dualProviralAlignedReads = defaultdict(list)
  hostReadsWithPotentialChimera = defaultdict(list)
  unmappedPotentialChimera = defaultdict(list)

  # recover all proviral "chromosome" names from partial fasta file used by Cellranger
  print("Getting proviral records")
  proviralSeqs = defaultdict(lambda: [])
  proviralFastaIds = getProviralFastaIDs(args.viralFasta, proviralSeqs)

  # get possible LTR regions from fasta file
  # TODO build into program (as opposed to running it before)
  print("Getting potential LTRs")
  potentialLTR = parseLTRMatches(args.LTRmatches, proviralSeqs)
  
  if not os.path.exists(args.outputDir + "/" + outputFNs["proviralReads"]):
    # parse BAM file
    print("Parsing cellranger BAM (namesorted)")
    parseCellrangerBam(bamfile = args.bamfile,
      proviralFastaIds = proviralFastaIds,
      proviralReads = dualProviralAlignedReads,
      hostReadsWithPotentialChimera = hostReadsWithPotentialChimera,
      unmappedPotentialChimera = unmappedPotentialChimera,
      top_n = args.topNReads) #debugging

    # output BAM files
    print("Writing out BAM files of parsed records")

    cellrangerBam = pysam.AlignmentFile(args.bamfile, "rb")
    writeBam(args.outputDir + "/" + outputFNs["proviralReads"], cellrangerBam, dualProviralAlignedReads)
    writeBam(args.outputDir + "/" + outputFNs["hostWithPotentialChimera"], cellrangerBam, hostReadsWithPotentialChimera)
    writeBam(args.outputDir + "/" +  outputFNs["umappedWithPotentialChimera"], cellrangerBam, unmappedPotentialChimera)
    cellrangerBam.close()

  else:
    print("Parsed BAM files already found. Importing these files to save time.")
    
    # import files...
    #dualProviralAlignedReads = importProcessedBam(args.outputDir + "/" + outputFNs["proviralReads"], returnDict = True)
    hostReadsWithPotentialChimera = importProcessedBam(args.outputDir + "/" + outputFNs["hostWithPotentialChimera"], returnDict = True)
    #unmappedPotentialChimera = importProcessedBam(args.outputDir + "/" +  outputFNs["umappedWithPotentialChimera"], returnDict = True)

  # parse host reads with potential chimera
  print("Finding valid chimeras from host reads")
  hostValidChimeras = parseHostReadsWithPotentialChimera(hostReadsWithPotentialChimera,
    potentialLTR,
    clipMinLen = args.LTRClipLen)
  
  cellrangerBam = pysam.AlignmentFile(args.bamfile, "rb")
  writeBam(args.outputDir + "/" + outputFNs["hostWithValidChimera"], cellrangerBam, hostValidChimeras["validReads"])
  print(hostValidChimeras["validHits"])
  cellrangerBam.close()


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

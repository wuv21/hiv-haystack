import pysam
from Bio import SeqIO
from collections import defaultdict
import logging
import argparse
import os


def getProviralFastaIDs(fafile):
  ids = []
  for record in SeqIO.parse(fafile, format = "fasta"):
    ids.append(record.id)

  return ids


def extractCellBarcode(read):
    tags = dict(read.tags)

    # accept only CB tag because it passes the allowlist set by 10X
    if "CB" in tags:
        barcode = tags["CB"]
    else:
        barcode = "None"

    return barcode


def parseCellrangerBam(bamfile, proviralFastaIds, proviralReads, hostReadsWithPotentialChimera, unmappedPotentialChimera, top_n = -1):
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 20)

  readIndex = 0
  for read in bam:
    refnameIsProviral = read.reference_name in proviralFastaIds
    # supposed to take mate's ref name or if no mate, the next record in BAM file
    nextRefnameIsProviral = read.next_reference_name in proviralFastaIds

    # ignore if optical/PCR duplicate OR without a mate
    if (read.flag & 1024) or (not read.flag & 1):
      readIndex += 1
      continue
    
    # if read is properly mapped in a pair AND not proviral aligned AND there is soft clipping involved
    # TODO add only if S number is > than settings
    elif (read.flag & 2) and (not refnameIsProviral) and ("S" in read.cigarstring):
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
  # initial fn arguments
  # TODO parameterize these into cmd line arguments
  
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
  proviralFastaIds = getProviralFastaIDs(args.viralFasta)

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
  args = parser.parse_args()

  if not os.path.exists(args.outputDir):
    os.makedirs(args.outputDir)

  main(args)

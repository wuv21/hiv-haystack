import pysam
from Bio import SeqIO
from collections import defaultdict

def getProviralFastaIDs(fafile):
  ids = []
  for record in SeqIO.parse(fafile):
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
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 10)

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
    elif (not read.flag & 14) and (read.reference_name in proviralFastaIds or read.next_reference_name in proviralFastaIds):
      # move to chimera identification
      unmappedPotentialChimera[read.qname].append(read)
    
    readIndex += 1
    if top_n != -1 and readIndex > top_n:
      return
    
  return bam



def main():
  fnArgs = {
    "cellranger_ns_bam": "hello",
    "bwa_proviral_fa": "hello",
    "bwa_host_fa": "hello"
  }

  dualProviralAlignedReads = defaultdict(list)
  hostReadsWithPotentialChimera = defaultdict(list)
  unmappedPotentialChimera = defaultdict(list)

  proviralFastaIds = getProviralFastaIDs(fnArgs["bwa_proviral_fa"])
  cellrangerBam = parseCellrangerBam(fnArgs["cellranger_ns_bam"],
    proviralFastaIds = proviralFastaIds,
    proviralReads = dualProviralAlignedReads,
    hostReadsWithPotentialChimera = hostReadsWithPotentialChimera,
    unmappedPotentialChimera = unmappedPotentialChimera,
    top_n = -1) #debugging


if __name__ == '__main__':
  main()
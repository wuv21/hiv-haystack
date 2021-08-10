import pysam
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
      val[read.query_name].append(read)
    else:
      val.append(read)

  return val


def writeFasta(chimeras, hostClipFastaFn):
  records = []
  for qnameKey in chimeras:
    chimera = chimeras[qnameKey]
    if chimera["adjustedHostSoftClip"] is not None:
      seq = Seq(chimera["adjustedHostSoftClip"])
    else:
      seq = Seq(chimera["hostSoftClip"]["clippedFrag"])
    
    record = SeqRecord(
      id = chimera["read"].qname,
      seq = seq,
      description = ""
    )

    records.append(record)

  SeqIO.write(records, hostClipFastaFn, "fasta")
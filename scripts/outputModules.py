from scripts.baseFunctions import extractCellBarcode, separateCigarString
from csv import writer

class IntegrationSite(object):
  def __init__(self, chr, orient, pos):
    super().__init__()

    self.chr = chr
    self.orient = orient
    self.pos = pos

  def __str__(self):
    return("int site at {}{}{}".format(self.chr, self.orient, self.pos))

  def returnAsList(self):
    return [self.chr, self.orient, self.pos]


class ProviralFragment(object):
  def __init__(self):
    super().__init__()
    
    self.seqname = ""
    self.startBp = 0 #0-based start pos
    self.endBp = 0 #0-based end pos (endBp is the actual end Bp as opposed to position + 1)
    self.cbc = ""
    self.readname = ""
    self.usingAlt = None
    self.confirmedAlt = False
    self.alreadyRecordedInIntegration = False

  def __str__(self):
    return "{} {}:{}-{}".format(self.cbc, self.seqname, self.startBp, self.endBp)

  def setManually(self, seqname, startBp, endBp, cbc, readname, usingAlt = None):
    self.seqname = seqname
    self.startBp = startBp #0-based start pos
    self.endBp = endBp #0-based end pos (endBp is the actual end Bp as opposed to position + 1)
    self.cbc = cbc
    self.usingAlt = usingAlt
    self.readname = readname

  def setFromRead(self, read):
    self.seqname = read.reference_name
    self.startBp = read.reference_start
    self.endBp = read.reference_end - 1
    self.cbc = extractCellBarcode(read)
    self.readname = read.qname

  def setIntegrationAnalysisFlag(self, status):
    self.alreadyRecordedInIntegration = status

  def setAlt(self, usingAlt):
    if usingAlt is None:
      self.usingAlt = None
    elif type(usingAlt) == list and len(usingAlt) == 0:
      self.usingAlt = None
    else:
      self.usingAlt = usingAlt

  def confirmedAltCase(self):
    newSeqname, newPos, newCigarstring = self.usingAlt
    # newCigar = separateCigarString(newCigarstring)
    origLen = self.endBp - self.startBp + 1

    self.seqname = newSeqname
    self.newPos = newPos - 1 # to follow 0-index since pysam doesn't auto adjust this, unlike for read.reference_start
    self.endPos = self.newPos + origLen - 1 # following same nomenclature

    self.confirmedAlt = True

  def returnAsList(self):
    return [self.cbc, self.seqname, self.startBp, self.endBp,
    self.readname, str(self.usingAlt), str(self.confirmedAlt), str(self.alreadyRecordedInIntegration)]

  
class ChimericRead(object):
  def __init__(self, read, intsite, proviralFragment):
    super().__init__()
    
    self.read = read
    self.intsite = intsite
    self.proviralFragment = proviralFragment

  def __str__(self):
    return "{} is chimeric with {}. Proviral fragment: {}".format(
      self.read.query_name,
      str(self.intsite),
      str(self.proviralFragment))


class ReadPairDualProviral(object):
  def __init__(self, read1 : ProviralFragment, read2 : ProviralFragment):
    super().__init__()
    self.read1 = read1
    self.read2 = read2
    self.potentialEditRead = ""
    self.potentialEditData = None
    self.potentialEditIsAlt = False

  def setPotentialClipEdit(self, readNum, readData, isAlt):
    self.potentialEditRead = readNum
    self.potentialEditData = readData
    self.potentialEditIsAlt = isAlt

  def unsetPotentialClipEdit(self):
    self.potentialEditRead = ""
    self.potentialEditData = None
    self.potentialEditIsAlt = False

  def updateWithConfirmedEdit(self, newProviralFrag):
    if self.potentialEditRead == "read1":
      self.read1 = newProviralFrag
      self.read2.confirmedAltCase()

    elif self.potentialEditRead == "read2":
      self.read2 = newProviralFrag
      self.read1.confirmedAltCase()

  def returnAsList(self):
    read1List = self.read1.returnAsList()
    read2List = self.read2.returnAsList()

    return [read1List, read2List]


class CompiledDataset(object):
  def __init__(self,
    validChimerasFromHostReads,
    validChimerasFromViralReads,
    validChimerasFromUnmappedReadsHost,
    validChimerasFromUnmappedReadsViral,
    validViralReads,
    unmappedViralReads):

    super().__init__()
    self.integrationSites = []
    self.pairedViralFrags = []
    self.collatedViralFrags = []

    if validChimerasFromViralReads is not None:
      for k in validChimerasFromViralReads:
        keypair = validChimerasFromViralReads[k]
        
        for c in keypair:
          self.integrationSites.append(c)
          self.validViralReads[c.proviralFragment.readname].setIntegrationAnalysisFlag(True)

          # self.collatedViralFrags.append(c.proviralFragment.returnAsList())

    for x in validChimerasFromHostReads:
      if len(x['minus']) != 0:
        self.integrationSites = self.integrationSites + x['minus']

      elif len(x['plus']) != 0:
        self.integrationSites = self.integrationSites + x['plus']

    for x in validChimerasFromUnmappedReadsHost:
      if len(x['minus']) != 0:
        self.integrationSites = self.integrationSites + x['minus']
        
        # TODO need to fix for multiple hits...
        self.collatedViralFrags.append(x['minus'][0].proviralFragment.returnAsList())

      elif len(x['plus']) != 0:
        self.integrationSites = self.integrationSites + x['plus']

        self.collatedViralFrags.append(x['plus'][0].proviralFragment.returnAsList())
    
    unmappedValidChimeraReadNames = []
    if validChimerasFromUnmappedReadsViral is not None:
      for key in validChimerasFromUnmappedReadsViral:
        alignedSites = validChimerasFromUnmappedReadsViral[key]
        for i in alignedSites:
          self.integrationSites.append(i)
          unmappedValidChimeraReadNames.append(i.proviralFragment.readname)

    # parse through paired viral reads
    for v in validViralReads:
      readPair = validViralReads[v]
      self.pairedViralFrags.append(readPair)
      
      readPairList = readPair.returnAsList()
      self.collatedViralFrags.append(readPairList[0])
      self.collatedViralFrags.append(readPairList[1])

    # parse through unampped viral reads
    for v in unmappedViralReads:
      if v.readname in unmappedValidChimeraReadNames:
        v.setIntegrationAnalysisFlag(True)
      
      self.collatedViralFrags.append(v.returnAsList())


  def exportIntegrationSiteTSV(self, fnIntSite, fnIntSiteFrag):
    output = [[x.proviralFragment.cbc] + x.intsite.returnAsList() for x in self.integrationSites]
    outputPV = [x.proviralFragment.returnAsList() for x in self.integrationSites]

    # export integration sites
    with open(fnIntSite, "w") as tsvfile:
      writ1 = writer(tsvfile, delimiter = "\t")

      writ1.writerow(["cbc", "chr", "orient", "pos"])
      for o in output:
        writ1.writerow(o)
    
    # export proviral frags from integration sites
    with open(fnIntSiteFrag, "w") as tsvfile2:
      writ2 = writer(tsvfile2, delimiter = "\t")

      writ2.writerow(["cbc", "seqname", "startBp", "endBp",
      "readname", "usingAlt", "confirmedAlt"])
      for o in outputPV:
        writ2.writerow(o[:-1])


  def exportProviralCoverageTSV(self, fn):
    with open(fn, "w") as tsvfile:
      writ = writer(tsvfile, delimiter = "\t")

      writ.writerow(["cbc", "seqname", "startBp", "endBp",
      "readname", "usingAlt", "confirmedAlt", "alreadyRecordedInIntegration"])
      for o in self.collatedViralFrags:
        writ.writerow(o)

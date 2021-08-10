from scripts.baseFunctions import extractCellBarcode

class IntegrationSite(object):
  def __init__(self, chr, orient, pos):
    super().__init__()

    self.chr = chr
    self.orient = orient
    self.pos = pos

  def __str__(self):
    return("int site at {}{}{}".format(self.chr, self.orient, self.pos))


class ProviralFragment(object):
  def __init__(self):
    super().__init__()
    
    self.seqname = ""
    self.startBp = 0 #0-based start pos
    self.endBp = 0 #0-based end pos (endBp is the actual end Bp as opposed to position + 1)
    self.cbc = ""
    self.usingAlt = None

  def __str__(self):
    return "{} {}:{}-{}".format(self.cbc, self.seqname, self.startBp, self.endBp)

  def setManually(self, seqname, startBp, endBp, cbc, usingAlt = None):
    self.seqname = seqname
    self.startBp = startBp #0-based start pos
    self.endBp = endBp #0-based end pos (endBp is the actual end Bp as opposed to position + 1)
    self.cbc = cbc
    self.usingAlt = usingAlt

  def setFromRead(self, read):
    self.seqname = read.reference_name
    self.startBp = read.reference_start
    self.endBp = read.reference_end - 1
    self.cbc = extractCellBarcode(read)

  def setAlt(self, usingAlt):
    if usingAlt is None:
      self.usingAlt = None
    elif type(usingAlt) == list and len(usingAlt) == 0:
      self.usingAlt = None
    else:
      self.usingAlt = usingAlt

  
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


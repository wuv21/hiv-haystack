class IntegrationSite(object):
  def __init__(self, chr, orient, pos):
    super().__init__()

    self.chr = chr
    self.orient = orient
    self.pos = pos

  def __str__(self):
    return("int site at {}{}{}".format(self.chr, self.orient, self.pos))


class ProviralFragment(object):
  def __init__(self, seqname, startBp, endBp, cbc, usingAlt = None):
    super().__init__()
    
    self.seqname = seqname
    self.startBp = startBp #0-based start pos
    self.endBp = endBp #0-based end pos (endBp is the actual end Bp as opposed to position + 1)
    self.cbc = cbc
    self.usingAlt = usingAlt

  def __str__(self):
    return "{} {}:{}-{}".format(self.cbc, self.seqname, self.startBp, self.endBp)

  
class ChimericRead(object):
  def __init__(self, read, intsite, proviralFragment):
    super().__init__()
    
    self.read = read
    self.intsite = intsite
    self.proviralFragment = proviralFragment

  def __str__(self):
    return "{} is chimeric read with {}. Proviral fragment: {}".format(self.read.qname, str(self.intsite), str(self.proviralFragment))


class ReadPairWithChimericRead(object):
  def __init__(self, chimericRead, nonChimericRead):
    super().__init__()
    self.chimericRead = chimericRead
    self.nonChimericRead = nonChimericRead

  def isNonChimericReadProviral(self):
    return type(self.nonChimericRead) is ProviralFragment 


class ReadPairDualProviral(object):
  def __init__(self, read1, read2):
    super().__init__()
    self.read1 = read1
    self.read2 = read2


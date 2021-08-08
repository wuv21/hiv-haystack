class IntegrationSite(object):
  def __init__(self, chr, pos):
    super().__init__()

    self.chr = chr
    self.pos = pos


class ProviralFragment(object):
  def __init__(self, seqname, startBp, endBp, orient, usingAlt = None):
    super().__init__()
    
    self.seqname = seqname
    self.startBp = startBp
    self.endBp = endBp
    self.orient = orient
    self.usingAlt = usingAlt

  
class ChimericRead(object):
  def __init__(self, read, intsite, proviralFragment):
    super().__init__()
    
    self.read = read
    self.intsite = intsite
    self.proviralFragment = proviralFragment


class ReadPairWithChimericRead(object):
  def __init__(self, chimericRead, nonChimericRead, cbc):
    super().__init__()
    self.chimericRead = chimericRead
    self.nonChimericRead = nonChimericRead
    self.cbc = cbc

  def isNonChimericReadProviral(self):
    return type(self.nonChimericRead) is ProviralFragment 


class ReadPairDualProviral(object):
  def __init__(self, read1, read2, cbc):
    super().__init__()
    self.read1 = read1
    self.read2 = read2
    self.cbc = cbc


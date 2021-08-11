import re

def separateCigarString(cigarstring):
  cigarSep = re.findall(r"(\d+\w)", cigarstring)
  cigarSepExpanded = [re.split(r"(\d+)", x)[1:3] for x in cigarSep]

  return cigarSepExpanded

def extractCellBarcode(read):
  # accept only CB tag because it passes the allowlist set by 10X
  tags = dict(read.tags)
  if "CB" in tags:
      barcode = tags["CB"]
  else:
      barcode = None

  return barcode

def getAltAlign(read):
  if not read.has_tag("XA"):
    return None

  altAlignRaw = read.get_tag("XA")
  
  # remove last semicolon
  altAlignRaw = altAlignRaw[:-1]
  altAligns = altAlignRaw.split(";")
  altAligns = [x.split(",") for x in altAligns]

  return altAligns

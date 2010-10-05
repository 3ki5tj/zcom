#!/usr/bin/env python
import os, sys, re 
from copy import copy

class CPreprocessor:
  '''
  preprocessor, taken literally, not understood
  properties:
    empty
    raw
    begin, end
    pif, pelse, pendif
  '''
  def __init__(self, src, p):
    self.empty = 1
    if src: self.parse(src, p)

  def __copy__(s):
    c = CPreprocessor(None, None)
    c.empty = s.empty
    c.raw = s.raw
    c.begin = copy(s.begin)
    c.end = copy(s.end)
    c.pif = s.pif
    c.pelse = s.pelse
    c.pendif = s.pendif
    c.flag = s.flag
    c.flagval = s.flagval
    c.cmt = s.cmt
    return c

  def __deepcopy__(s, memo):
    return copy(s)

  def parse(self, src, p):
    s = p.skipspace(src)
    self.begin = copy(p)

    if not s.startswith("#"): return -1
    self.raw = s[1:].strip()
    self.empty = 0
    p.col += len(s)
    self.end = copy(p)

    self.pif = None
    self.pelse = 0
    self.pendif = 0
    self.flag = self.flagval = self.cmt = None

    # try to understand what the preprocessor does
    s = self.raw
    if s.startswith("endif"):
      self.pendif = 1
    elif s.startswith("if"):
      if s.startswith("ifdef"):
        self.pif = "defined(%s)" % s[5:].strip()
      elif s.startswith("ifndef"):
        self.pif = "!defined(%s)" % s[6:].strip()
      else:
        self.pif = s[2:].strip()
    elif s.startswith("else"):
      self.pelse = 1
    elif s.startswith("elif"):
      self.pelse = 1
      self.pif = s[4:].strip()
    elif s.startswith("define"):
      self.doflag(src, p)
    else:
      pass

  def doflag(self, src, p):
    s = self.raw
    pat = r"define\s+([a-zA-Z_]\w*)\s+(\w*)\s*"
    m = re.match(pat, s)
    if not m: return
    self.flag = m.group(1)
    self.flagval = m.group(2)
    rest0 = s[m.end(2):]
    rest = rest0.strip()
    if rest.startswith("//"):
      self.cmt = rest
    elif rest.startswith("/*") and rest.endswith("*/"):
      self.cmt = rest

    if self.cmt:
      # fork over the comment part
      p.col -= len(rest0)
      self.end = copy(p)
      self.raw = s[:m.end(0)]
      #print "%s: %s" % (p, p.getline(src)); raw_input()

  def isempty(self): return self.empty



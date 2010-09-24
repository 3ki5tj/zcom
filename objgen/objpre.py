#!/usr/bin/env python
import os, sys, re 
from copy import copy

class CPreprocessor:
  '''
  preprocessor, taken literally, not understood
  '''
  def __init__(self, src, p):
    self.empty = 1
    self.parse(src, p)

  def parse(self, src, p):
    s = p.skipspace(src)
    self.begin = copy(p)

    if not s.startswith("#"): return -1
    self.raw = s[1:].strip()
    self.empty = 0
    p.col += len(s)

    self.pif = None
    self.pelse = 0
    self.pendif = 0

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
    else:
      pass

  def isempty(self): return self.empty



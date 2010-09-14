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

  def isempty(self): return self.empty



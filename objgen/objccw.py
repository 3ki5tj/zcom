#!/usr/bin/env python
import os, sys, re
from copy import copy

class CCodeWriter:
  '''
  output code in an indented fashion
  '''
  sindent = "  "
  def __init__(self, nindents = 0):
    self.nindents = nindents
    self.s = ""

  def inc(self): self.nindents += 1
  def dec(self): 
    self.nindents = self.nindents-1 if self.nindents > 0 else 0

  def add(self, t, *args):
    t = t % args
    if t.lstrip().startswith("}"): self.dec()
    if self.s.endswith("\n") and not t.startswith("#"):
      self.s += self.sindent * self.nindents
    self.s += t
    if t.rstrip().endswith("{"): self.inc()

  def addln(self, t, *args):
    self.add( (t.strip() % args) + '\n' )

  def remove_empty_pp(self):
    ''' remove empty preprocessor blocks 
    #if XXX
    #else
    #endif
    '''
    lines = self.s.splitlines()
    i = 1
    while i < len(lines):
      if lines[i].startswith("#endif") and i > 0:
        lnprev = lines[i-1]
        if (lnprev.startswith("#else") or 
            lnprev.startswith("#elif")):
          # remove an empty #else/#elif branch
          lines = lines[:i-1] + lines[i:]
          continue
        elif lnprev.startswith("#if"):
          lines = lines[:i-1] + lines[i+1:]
          continue
      i += 1
    self.s = '\n'.join(lines)

  def gets(self):
    self.remove_empty_pp()
    return self.s
    



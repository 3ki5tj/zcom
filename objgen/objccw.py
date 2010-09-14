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

  def gets(self): return self.s
    



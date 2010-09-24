#!/usr/bin/env python
'''
C comment

TODO:
  * comment fusion
'''

import os, sys, re 
from copy import copy

cmt0 = "/*"
cmt1 = "*/"
lncmt = "//"

class CComment:
  '''
  C block-comment or one-line comment
  '''
  def __init__(self, src, p, maxnl = 0):
    self.empty = 1
    if not src or 0 != self.span(src, p, maxnl):    # .begin, .end, and .raw
      return
    self.empty = 0

  def __copy__(s):
    c = CComment(None, None)
    c.empty = s.empty
    c.raw = s.raw
    c.begin = copy(s.begin)
    c.end = copy(s.end)
    return c

  def __deepcopy__(self, memo):
    return copy(self)
    
  def span(self, src, p, maxnl = 0):
    '''
    find the beginning and end of a comment,
    starting from pos
    '''
    s = p.skipspace(src, maxnl) # note, its skips blank lines if maxnl = 0
    if s == None: return -1

    if s.startswith(cmt0) or s.startswith(lncmt):
      self.type = "block" if s.startswith(cmt0) else "line"
      self.begin = copy(p)
    else:
      #print "missing comment beginning mark %s" % p
      return -1

    p.col += 2 # skip `/*' or `//'

    if s.startswith(lncmt): # line comment
      self.raw = s[2:].strip()
      p.endofline(src)
      self.end = copy(p)
    else: # block comment
      # aggressively search the end of comment
      while 1:
        s = p.skipspace(src)
        if s.startswith(cmt1):
          p.col += len(cmt1)
          self.end = copy(p)
          #print "comment found from %s to %s" % (self.begin, self.end)
          break
        p.gettoken(src) # get token by token
      else:
        print "Warning: cannot find the end of comment from %s, now = %s" % (
            self.begin, now)
        raise Exception
      self.raw = self.extract_text(src)
    return 0

  def extract_text(self, src):
    '''
    extract text within block comment
    '''
    self.checkbe()
    s = self.begin.gettext(src, self.end).strip()
    if s.startswith(cmt0) and s.endswith(cmt1):
      s = s[len(cmt0) : len(s) - len(cmt1)]
    #print "comment text is [%s] from %s to %s" % (s, self.begin, self.end)
    return s

  def checkbe(self):
    self.begin.check()
    self.end.check()

  def isempty(self): return self.empty
  def __str__(self): return self.raw
  def __repr__(self): return __str__(self)




#!/usr/bin/env python
import os, sys, re
from copy import copy

'''
simple objects 
'''

class P:
  '''
  position in the code (line, column)
  unfortunately, the most vulnerable class
  '''
  def __init__(self, row = 0, col = 0):
    self.row = row
    self.col = col
    self.check() # see if input are integers
  def __str__(self):
    self.check()
    return "line: %d, col %d" % (self.row + 1, self.col + 1)
  def __repr__(self):
    self.check();
    return "%d, %d" % (self.row, self.col)
  def __cmp__(a, b):
    a.check()
    if type(a) != type(b):
      print "trying to compare %s with %s" % (type(a), type(b))
      raise Exception
    b.check()
    return (a.row - b.row) * 10000 + a.col - b.col
  def check(self):
    if type(self.row) != int:
      print "row %s is not a number" % (self.row)
      raise Exception
    if type(self.col) != int:
      print "col %s is not a number" % (self.col)
      raise Exception
  def getline(p, src):
    ''' get a line starting from the current position '''
    p.check()
    return src[p.row][p.col:] if p.row < len(src) else None
  def endofline(p, src):
    ''' move to the end of this line '''
    p.check()
    s = p.getline(src)
    if s: p.col += len(s)
  def nextline(self, src = None):
    '''
    go to the next line
    we don't really need src, just make it look like other methods
    '''
    self.check()
    self.row += 1
    self.col = 0
  def skipspace(p, src, maxnl = 0):
    '''
    skip spaces (including multiple blank lines)
    return the line
    `maxnl' specifies the maximal number of newlines to skip
      0: no limit, 
      1: to the end of this line, 
      2: this line and the next line
      ...
    '''
    p.check()
    i = 0
    while 1:
      s = p.getline(src)
      if s == None: # end of file
        self.ttype = self.token = None
        break
      # skip leading space
      m = re.match(r"(\s*).*", s)
      if not m: raise Exception
      p.col += m.end(1)
      s = p.getline(src)
      if len(s) != 0: break # found a token
      i += 1
      if maxnl > 0 and i > maxnl: break
      # this line is exhausted
      p.nextline()
      continue
    return s

  def gettoken(p, src, maxnl = 0):
    '''
    aggressively get a token starting from the current `p'
    `p' is updated to the end of the current token
    return the token only
    '''
    s = p.skipspace(src, maxnl)
    
    # try to get a word
    if s != None:
      m = re.match(r"([a-zA-Z_]\w*).*", s)
      if m:
        p.col += m.end(1)
        p.ttype, p.token = "word", s[:m.end(1)]
      else:
        p.col += 1
        p.ttype = p.token = s[0]
    return p.token, p.ttype

  def ungettok(p, src):
    '''
    try to unget a token, it will not unget space
    '''
    p.check()
    tok = p.token
    if not tok:
      print "cannot unget None"
      raise Exception
    p.col -= len(tok)
    if (p.col < 0 or
        not p.getline(src).startswith(tok)):
      raise Exception
    p.token = p.ttype  = None

  def peektok(p, src):
    '''
    peek at the next token
    '''
    p.check()
    q = copy(p) # create a backup
    return q.gettoken(src)

  def gettext(self, src, q):
    '''
    get entire string from p to q
    '''
    p = copy(self) # so self is intact
    if q <= P(0,0):
      q = P(len(src), 0) # end of string
    s = ""
    while p.row < q.row:
      s += p.getline(src)
      p.nextline()
    if p.row == q.row and p.col < q.col:
      s += src[p.row][p.col: q.col]
    return s




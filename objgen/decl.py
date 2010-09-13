#!/usr/bin/env python
import os, sys, re, filecmp

'''
python implementation of K & R's  dcl program
'''

class CDecl:
  def __init__(self, s):
    if type(s) != type(""): raise Exception
    self.raw = s.strip()
    self.nraw = len(self.raw)
    self.types = []
    self.dcl_wrapper()

  def __str__(self):
    return ' '.join(self.types)

  def dcl_wrapper(self):
    '''
    get the data type, then call dcl
    '''
    self.pos = 0
    self.token = self.ttype = None
    self.gettoken()
    if self.ttype != "word":
      print "complex type is yet supported, [%s]" % self.raw
      print self.ttype, self.token, self.pos
      raise Exception
    self.datatype = self.token
    self.dcl()
    self.types += [self.datatype]

  def dcl(self):
    '''
    declaration
    '''
    ns = 0
    while self.pos < self.nraw:
      self.gettoken()
      if self.ttype != "*": break
      ns += 1
    self.dirdcl()
    self.types += ["pointer to"] * ns

  def dirdcl(self):
    '''
    direct declaration
    '''
    if self.ttype == '(': # ( dcl )
      self.dcl()
      if self.ttype != ')':
        print "missing ), s: %s, pos %d" % (self.raw, self.pos)
        raise Exception
    elif self.ttype == "word":
      self.name = self.token
    else:
      print "expected name or ( at pos: %d, s: %s" % (self.pos, self.raw)
    while 1:
      self.gettoken()
      if self.ttype == "()":
        self.types += ["function returning"]
      elif self.ttype == "[]":
        self.types += ["array " + self.token + " of"]
      else: break

  def gettoken(self):
    '''
    get a token starting from the current pos
    '''
    s = self.raw[self.pos:]
    m = re.match(r"(\s*).*", s)
    if not m: raise Exception
    self.pos += m.end(1)
    s = self.raw[self.pos:]

    if len(s) == 0:
      self.ttype = self.token = None
      return
    m = re.match(r"([a-zA-Z_]\w*).*", s)
    if m:
      self.ttype, self.token = "word", s[:m.end(1)]
      self.pos += m.end(1)
      return
    m = re.match(r"(\[)(.*?)(\]).*", s)
    if m:
      self.ttype, self.token = "[]", s[m.start(2):m.end(2)]
      self.pos += m.end(3)
      return
    if s.startswith("()"):
      self.ttype = self.token = "()"
      self.pos += 2;
      return
    self.ttype = self.token = s[0]
    self.pos += 1

def main():
  print CDecl(sys.argv[1])
 
if __name__ == "__main__":
  main()


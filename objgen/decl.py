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
      if self.gettoken() != "*": break
      ns += 1
    self.dirdcl()
    self.types += ["pointer to"] * ns

  def dirdcl(self):
    '''
    direct declaration
    '''
    if self.ttype == '(': # ( dcl )
      self.dcl() # in simplest case, dcl() reads a word, e.g. "abc",
                 # then calls dirdcl() again, which copies the word 
                 # to self.name, see below, and reads another token
                 # in the while 1 loop, in which ")" should be read
                 # and returns to dcl(). dcl() prints the number of
                 # stars and returns to here, so the last token
                 # should still be ")" after two levels of recursion
                 # the recursion is necessary because, dcl() really
                 # does nothing but count the number of stars
      if self.ttype != ')':
        print "missing ), s: %s, pos %d" % (self.raw, self.pos)
        raise Exception
    elif self.ttype == "word":
      self.name = self.token
    else:
      print "expected name or ( at pos: %d, s: %s" % (self.pos, self.raw)
    while 1:
      ttype = self.gettoken()
      if ttype == "(":
        self.ptlist()
        if self.ttype != ')':
          print "missing ) after parameter type list"
        self.types += ["function returning"]
      elif ttype == "[]":
        self.types += ["array " + self.token + " of"]
      else: 
        break 

  def ptlist(self):
    self.gettoken()
    #if self.ttype != ')':
    #  print "cannot handle complex function"
    #  raise Exception

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
      return self.ttype
    m = re.match(r"([a-zA-Z_]\w*).*", s) # a variable
    if m:
      self.ttype, self.token = "word", s[:m.end(1)]
      self.pos += m.end(1)
      return self.ttype
    m = re.match(r"(\[)(.*?)(\]).*", s)  # [...]
    if m:
      self.ttype, self.token = "[]", s[m.start(2):m.end(2)]
      self.pos += m.end(3)
      return self.ttype
    self.ttype = self.token = s[0]
    self.pos += 1
    return self.ttype

def main():
  print CDecl(sys.argv[1])
 
if __name__ == "__main__":
  main()


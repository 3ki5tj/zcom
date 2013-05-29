#!/usr/bin/env python

''' routines to handle C source code
    Cheng Zhang (c) 2013

  class CParser
    * switchpp():       switch the current preprocessor state
    * switchcmt():      switch the current comment state
    
  Standalone helper function
    * endingcmt():      if a code line starts a block comment

'''


import re


def endingcmt(ln, incmt = False):
  ''' check if this line ends with a block comment
        a=b; /* hahaha */ b=3; /* comment
      NOTE: it will fail if `/*' or '*/' is hidden in a string
  '''

  # finish the current comment
  if incmt:
    # currently in a block comment
    m = re.match(r".*?\*/(.*)", ln)
    if m:
      ln = m.group(1)
    else: # the block comment is not done
      return True
 
  #print "input:", ln 
  while 1:
    m = re.match(r".*?(/\*.*?\*/)(.*)", ln)
    m1 = re.match(r".*?(/\*).*", ln)
    m2 = re.match(r".*?(//).*", ln)
    
    # check `//'
    if m2:
      if not m1: return False

      # both '//' and '/*' exist
      if m2.start(1) < m1.start(1):
        # `//' preceeds `/*'
        #print "line comment:", ln[m2.start(1):]
        return False
      elif not m:
        # `/*' preceeds `//', but no '*/'
        return True
      # closed `/* ... */', pass through
 
    if not m:
      # no more completed `/* ... */'
      #print "ending:", ln
      return (m1 != None)
    else:
      ln = m.group(2)


#print endingcmt("/*ha//h*/a/*//*/HA/*/*/HAHA")
#exit(1)


class CParser:
  ppn = 0  # preprocessor level
  inpp = False  # currently in a preprocessor
  incmt = False

  def __init__(self, ppn = 0, inpp = False, incmt = False):
    self.ppn = ppn
    self.inpp = inpp
    self.incmt = incmt


  def switchpp(self, ln):
    ''' switch the current state for preprocessors 
        according to the input line `ln'
        return if the current line is in a preprocessor block '''

    ln = ln.strip()
    # ignore the preprocessor, if it is in a comment
    if not self.incmt:
      if ln.startswith("#"):
        self.inpp = True
        if ln.startswith("#if"):
          self.ppn += 1
        elif ln.startswith("#endif"):
          self.ppn = max(self.ppn - 1, 0)

    inppthis = self.inpp

    # compute the pp state of the next line
    if self.inpp:
      # end the processor, unless there's a line continuation
      if not ln.endswith("\\"):
        self.inpp = False
    return inppthis


  def switchcmt(self, ln):
    ''' switch the current state for block comments
        according to the input line `ln'
        an embedded comment does not counts:
          a = 3; /* assignment */ b = 4;
        return if the current line is in a comment block '''

    ln = ln.strip()
    incmtthis = self.incmt
    if not self.incmt and ln.startswith("/*"):
      incmtthis = True
    self.incmt = endingcmt(ln, self.incmt)
    return incmtthis



def printppcmt(fn):
  ''' test: print preprocessor or comment lines '''
  s0 = open(fn).readlines()
  cp = CParser()
  for i in range(len(s0)):
    inpp = cp.switchpp(s0[i])
    incmt = cp.switchcmt(s0[i])
    if inpp or incmt:
      print "%5d:%d,%5s|%5s|%s" % (i+1, cp.ppn, cp.inpp, cp.incmt, s0[i]),



if __name__ == "__main__":
  import sys
  if len(sys.argv) > 1: printppcmt( sys.argv[1] )



#!/usr/bin/env python

''' routines to handle C source code
    Cheng Zhang (c) 2013

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



def switchpp(ln, pplevel, inpp, incmt):
  ''' switch the current state for preprocessors '''

  ln = ln.strip()

  # ignore the preprocessor, if it is in a comment
  if not incmt:
    if ln.startswith("#"):
      inpp = True
      if ln.startswith("#if"):
        pplevel += 1
      elif ln.startswith("#endif"):
        pplevel = max(pplevel - 1, 0)

  inppthis = inpp

  # compute the pp state of the next line
  if inpp:
    # end the processor, unless there's a line continuation
    if not ln.endswith("\\"):
      inpp = False

  return pplevel, inpp, inppthis



def switchcmt(ln, incmt):
  ''' switch the current state for block comments 
      an embedded comment does not counts:
        a = 3; /* assignment */ b = 4;    
  '''

  ln = ln.strip()
  incmtthis = incmt
  if not incmt and ln.startswith("/*"):
    incmtthis = True
  
  return endingcmt(ln, incmt), incmtthis






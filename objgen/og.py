#!/usr/bin/env python

import os, sys, re, filecmp, glob, tempfile, shutil
from objgen import Parser, manual


def handle(fnout, template = None):
  ''' generate a C file according to a template
  the file name `template', if not given, is deduced from `fnout'
  e.g., fnout: `a.c',  template: `a.0.c' '''

  # guess the template name if it's empty
  if template == None:
    pt = os.path.splitext(fnout)
    template = pt[0] + '.0' + pt[1]
    ref = pt[0] + '.1' + pt[1]

  tmpfile = tempfile.mkstemp(suffix = pt[1], dir = ".")[1]

  # read in the template
  src = open(template, 'r').readlines()
  # construct a parser
  psr = Parser(src)
  # write the output
  open(tmpfile, 'w').write(psr.output())

  # if there is a reference file for the output
  if os.path.exists(ref):
    if not filecmp.cmp(tmpfile, ref):
      print "%s differs from %s" % (fnout, ref)

  # if `fnout' already exists, we back it file with a suffix
  #   .bak, .bak.1, .bak.2, ...
  if os.path.exists(fnout):
    if not filecmp.cmp(tmpfile, fnout): # different files
      fn = fnout + ".bak"
      i = 1
      while os.path.exists(fn):
        fn = fnout + ".bak." + str(i)
        i += 1

      shutil.copy2(fnout, fn)
    else:
      print "keep %s" % fnout
      os.remove(tmpfile)
      return

  # write the output file
  print "update %s" % fnout
  shutil.copy(tmpfile, fnout)
  os.remove(tmpfile)

def main():
  if "--help" in sys.argv or "-h" in sys.argv:
    print manual
    exit(0)

  if len(sys.argv) > 1:
    files = [sys.argv[1]]
  else:
    files = [re.sub(r"\.0", "", file) for file in glob.glob("*.0.[ch]")]
    #raw_input("about to handle: " +str(files)+ ", okay?")

  for fn in files:
    handle(fn)

if __name__ == "__main__":
  try:
    import psyco
    psyco.full()
  except ImportError: pass
  main()


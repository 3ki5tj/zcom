#!/usr/bin/env python

import os, sys, re, filecmp, glob, tempfile, shutil
from objgen import Parser


def handle(file, template = ""):
  ''' generate a C file according to a template '''
  # guess the template name if it's empty
  if template == "":
    pt = os.path.splitext(file)
    template = pt[0] + '.0' + pt[1]
    ref = pt[0] + '.1' + pt[1]
  
  tmpfile = tempfile.mkstemp(suffix = pt[1], dir = ".")[1]
  #raw_input(tmpfile)

  # read in the template
  src = open(template, 'r').readlines()
  # construct a parser
  psr = Parser(src)
  # write the processed file
  open(tmpfile, 'w').write(psr.output())

  if os.path.exists(ref):
    if not filecmp.cmp(tmpfile, ref):
      print "File %s is different from the reference %s" % (file, ref)

  if os.path.exists(file):
    if not filecmp.cmp(tmpfile, file): # different files
      # make a backup of the original file
      fn = file + ".bak"
      i = 1
      while os.path.exists(fn):
        fn = file + ".bak." + str(i)
        i += 1
      shutil.copy2(file, fn)
    else:
      print "no change is needed to %s" % file
      os.remove(tmpfile)
      return

  # write the output file
  print "update %s" % file
  shutil.copy(tmpfile, file)
  os.remove(tmpfile)

def main():
  if len(sys.argv) > 1: 
    files = [sys.argv[1]]
  else:
    files = [re.sub(r"\.0", "", file) for file in glob.glob("*.0.[ch]")]
    #raw_input("about to handle: " +str(files)+ ", okay?")

  for file in files:
    handle(file)

if __name__ == "__main__":
  try:
    import psyco
    psyco.full()
  except ImportError: pass
  main()


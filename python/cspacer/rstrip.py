#!/usr/bin/env python

''' remove trailing spaces of a text file
    Copyright (c) 2013 Cheng Zhang '''

import os, sys, getopt

verbose = 0
bytessaved = 0



def rtrimf(fn):
  ''' remove trailing spaces of file `fn' '''

  try:
    s0 = open(fn).read()
  except Exception:
    print "cannot open %s" % fn
    return

  # determine the end of line
  endl = '\n'
  lines = s0.splitlines()
  if len(lines):
    ln = lines[0]
    if ln.endswith("\r\n"):
      endl = "\r\n"
    elif ln.endswith("\n\r"):
      endl = "\n\r"
    else:
      endl = "\n"

  # remove trailing spaces
  s1 = [ln.rstrip() + endl for ln in lines]
  s1 = ''.join( s1 )

  bsaved = len(s0) - len(s1)
  global bytessaved
  bytessaved += bsaved

  if bsaved <= 0:
    if verbose >= 2:
      print "keeping %s" % fn
  else:
    if verbose >= 1:
      print "overwriting %s, len %s -> %s (save %s bytes)" % (
        fn, len(s0), len(s1), bsaved)
    open(fn, "w").writelines(s1)




def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] input"
  print """
  Remove trailing spaces

  OPTIONS:

   -l: include symbolic links
   -R: recursively apply to subdirectories, if `input' is a wildcard pattern
       like *.c, the pattern must be quoted as '*.c'
   -v: be verbose
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "lRhv",
         ["links", "recursive", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose

  recur = links = False
  for o, a in opts:
    if o in ("-R", "--recursive",):
      recur = True
    elif o in ("-l", "--links",):
      links = True
    elif o in ("-v",):
      verbose = 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  ls = args
  if not ls:
    # common text files
    pats = """*.c *.cpp *.h *.hpp *.java
              *.py *.pl *.rb *.php *.js
              *.f *.f90 *.f95 *.pas *.bas
              *.m *.ma *.gp *.tcl
              *.txt *.tex *.html *.htm
              *.cfg *.mdp
              *.sh *.csh
              *.bib
              *.md
              README* *akefile"""

    try: # limit the dependence on zcom.py
      import zcom
      ls = zcom.argsglob(args, pats, recur = recur, links = links)
    except ImportError:
      import glob
      ls = []
      for pat in pats.split():
        ls += glob.glob(pat)
  return ls



if __name__ == "__main__":
  fns = doargs()
  for fn in fns:
    rtrimf(fn)
  print "saved %s bytes" % bytessaved


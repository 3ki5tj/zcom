#!/usr/bin/env python

import os, sys, glob, getopt

def trim(fn):
  ''' remove trailing spaces of a file '''
  try:
    s0 = open(fn).read()
  except Exception:
    print "cannot open %s" % fn
    return
  s1 = ''.join( [ln.rstrip() + '\n' for ln in s0.splitlines(True)] )
  if s0 == s1:
    print "keep %s" % fn
  else:
    print "overwrite %s, len %s -> %s (save %s bytes)" % (
        fn, len(s0), len(s1), len(s0) - len(s1))
    open(fn, "w").writelines(s1)


def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] input_files"
  print "OPTIONS:"
  print " -R: apply to subdirectories"
  exit(1)


def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "Rh",
         ["recursive", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  recur = False
  for o, a in opts:
    if o in ("-R", "--recursive",):
      recur = True
    elif o in ("-h", "--help",):
      usage()

  ls = []
  if not recur: # list files under this dir
    for a in args:
      ls += glob.glob(a)
  else: # recursively list all files under all subdirectories
    root = os.getcwd()
    for r, d, f in os.walk("."):
      os.chdir(r)
      ls0 = []
      for a in args:
        ls0 += glob.glob(a)
      for fn0 in ls0:
        ls += [ os.path.join(r, fn0), ]
      os.chdir(root)

  if len(ls) <= 0:
    print "no files"
  return ls

def main():
  fns = doargs()
  for fn in fns:
    trim(fn)


if __name__ == "__main__":
  main()

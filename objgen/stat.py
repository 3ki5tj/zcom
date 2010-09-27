#!/usr/bin/env python
import os, sys, glob

def cntsize(list):
  nfiles = len(list)
  if not nfiles: return
  nlines = 0
  nchars = 0
  maxwid = max(len(f) for f in list)
  for f in list:
    nc = len(open(f).read())
    nl = len(open(f).readlines())
    print "%-*s: %6d lines, %8d chars, %6.2f chars per line" % (maxwid, f, nl, nc, 1.0 *nc/nl)
    nchars += nc
    nlines += nl
  print "%d files, %d lines, %d characters" % (nfiles, nlines, nchars)
  print "%.2f lines per file, %s characters per line" % (
      1.0 * nlines / nfiles, 1.0 * nchars / nlines)


if __name__ == "__main__":
  list = glob.glob("obj*.py") + ["og.py"]
  cntsize(list)
  


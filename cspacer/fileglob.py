#!/usr/bin/env python

''' recursively list all files matching a list of patterns
    Copyright (c) 2013 Cheng Zhang '''

import os, sys, glob

# skip hidden files
skiphidden = True



def fglob(pat, links = False, dir = None):
  ''' list files that match `pat' under `dir' '''

  if dir: pat = os.path.join(dir, pat)
  ls = [ a for a in glob.glob(pat) if os.path.isfile(a) ]
  if not links: # exclude symbolic links
    ls = [ a for a in ls if not os.path.islink(a) ]
  return ls



def fileglob(pats, links = False, recur = False):
  ''' find all files that match a list of patterns `pats' '''

  ls = []
  if not recur: # list files under this dir
    ls += [ g for pat in pats for g in fglob(pat, links) ]
  else: # recursively list all files under all subdirectories
    root = os.getcwd()
    for r, d, f in os.walk(root):

      if skiphidden:
        # stay away from hidden directories like `.git'
        danger = [ len(p) >= 2 and p[0] == "."
                 for p in r.split(os.sep) ]
        # treat booleans as 0 or 1, and sum them up
        if sum(danger) > 0: continue

      ls0 = [ g for p in pats for g in fglob(p, links, r) ]
      for fn0 in ls0:
        ls += [ os.path.join(r, fn0), ]

  # remove dupliated files by convert the list to a set
  # then convert it back to a list
  ls = list( set(ls) )
  return ls



def globargs(args, defpat = "*", links = False, recur = False):
  ''' glob patterns in the command line arguments arguments '''

  pats = defpat.split()
  if len(args) > 0:
    # parse the pattern in each argument
    pats = [ a for pat in args for a in pat.split() ]
  return fileglob(pats, True, recur)


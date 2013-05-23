#!/usr/bin/env python

import os, sys, glob, getopt

verbose = 0
bytessaved = 0

def trim(fn):
  ''' remove trailing spaces of file `fn' '''

  try:
    s0 = open(fn).read()
  except Exception:
    print "cannot open %s" % fn
    return

  # remove trailing spaces
  s1 = [ln.rstrip() + '\n' for ln in s0.splitlines()]
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
  OPTIONS:

   -s: include symbolic links
   -R: apply to subdirectories, if `input' is a wildcard pattern
       like *.c, the pattern must be quoted as '*.c'
   -v: be verbose
  """
  exit(1)


def fglob(pat, links = False, dir = None):
  ''' list files that match `pat' under `dir' '''

  if dir: pat = os.path.join(dir, pat)
  ls = [ a for a in glob.glob(pat) if os.path.isfile(a) ]
  if not links: # exclude symbolic links
    ls = [ a for a in ls if not os.path.islink(a) ]
  if verbose >= 3 and len(ls): print ls, pat
  return ls


def fileglob(pats, links = False, recur = False):
  ''' find all files that match a list of patterns `pats' '''

  ls = []
  if not recur: # list files under this dir
    ls += [ g for pat in pats for g in fglob(pat, links) ]
  else: # recursively list all files under all subdirectories
    root = os.getcwd()
    for r, d, f in os.walk(root):
      # stay away from `.git' ...
      danger = [ len(p) >= 2 and p[0] == "."
                 for p in r.split(os.sep) ]
      if sum(danger) > 0: continue
      ls0 = [ g for p in pats for g in fglob(p, links, r) ]
      for fn0 in ls0:
        ls += [ os.path.join(r, fn0), ]

  # remove dupliated files by convert the list to a set
  # then convert it back to a list
  ls = list( set(ls) )
  return ls



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

  # common text files
  pats = '''*.c *.cpp *.h *.hpp *.java *.py
            *.m *.ma *.gp *.txt *.tex
            *.html *.htm *.cfg *.mdp README* *akefile'''.split()
  if len(args) > 0:
    # parse the pattern in each argument
    pats = [ a for pat in args for a in pat.split() ]

  # compile a list of files
  ls = fileglob(pats, links, recur)
  if len(ls) <= 0: print "no file"
  return ls



def main():
  fns = doargs()
  for fn in fns: trim(fn)
  print "saved %s bytes" % bytessaved


if __name__ == "__main__":
  main()

#!/usr/bin/env python

''' remove unused functions from a python script
    Copyright (c) 2013 Cheng Zhang '''

import os, sys, re, getopt



verbose = 0
doall = False



def getfuncls(s):
  ''' collect functions '''

  ls = []
  for i in range(len(s)):
    ln = s[i].strip()
    m = re.match(r"def\s+(\w+?)\s*?\(.*", ln)
    if m: ls += [ m.group(1), ]
  return ls



def rmfunc(src, func, force):
  ''' remove unused functions from fn '''

  # 1. find the function definition block
  i0 = i1 = -1
  for i0 in range(len(src)):
    if re.search(r"^def\s+%s\s*\(" % func, src[i0]):
      break
  else:
    if verbose:
      print "cannot find the function block %s" % func
    return src

  for i1 in range(i0 + 1, len(src)):
    ln = src[i1]
    if ln.strip() != "" and not ln[0].isspace():
      break
  # skip blank lines
  while i1 < len(src) and src[i1].strip() == "":
    i1 += 1

  if verbose >= 2:
    print "function block for %s is %s" % (func, [i0, i1])

  # 2. remove the function
  rm = True
  if not force:
    # search occurrence
    for i in range(len(src)):
      if i >= i0 and i < i1: continue
      # if the function is used, we don't remove it
      if func + "(" in src[i]:
        rm = False
        break

  if rm:
    print "removing function '%s'" % func
    if verbose:
      print "function block for %s is\n%s" % (func, ''.join(src[i0:i1]) )

  # remove it
  if rm:
    src = src[:i0] + src[i1:]
  return src



def rmfuncs(src, funcls = None, force = False):
  if funcls == None: funcls = getfuncls(src)

  # iterate several times, for some functions become
  # unused only after other functions are removed
  niter = 1
  if not force: niter = len(funcls) + 1
  for i in range(niter):
    if verbose: print "rmfuncs: starting iteration %d" % i
    src0 = src[:]
    for func in funcls:
      src = rmfunc(src, func, force)
    if ''.join(src0) == ''.join(src):
      if verbose: print "rmfuncs: finishing at iteration %s" % i
      break
  return src



def usage():
  print """  Usage:
    %s [OPTIONS] input.py
  Remove unused functions of `input.py'
  Copyright (c) 2013 Cheng Zhang

  OPTIONS:
  -f, --func=       list all functions to be removed, separated by `,'
  -a, --all         remove all unused functions
  -o, --output=     output file
  -F, --force       remove even used functions
  -v                increase the verbosity
  --verbose=        set the verbosity
  -h, --help        help
  """ % sys.argv[0]
  exit(1)



def doargs():
  ''' handle the command-line arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "af:o:Fvh",
         ["all", "func=", "output=", "force", "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err)
    usage()

  global doall, verbose

  verbose = 0
  force = False
  fninp = fnout = None
  funcls = []

  for o, a in opts:
    if o in ("-a", "--all",):
      doall = True
    elif o in ("-f", "--func",):
      funcls = [ s.strip() for s in a.strip().split(",") ]
    elif o in ("-o", "--output",):
      fnout = a
    elif o in ("-F", "--force",):
      force = True
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  if len(args) > 0:
    fninp = args[0]
  else:
    print "need an input file"
    usage()

  if not os.path.exists(fninp):
    print "%s does not exist" % fninp
    usage()

  if not funcls and doall:
    funcls = None

  return fninp, fnout, funcls, force



if __name__ == "__main__":
  fninp, fnout, funcls, force = doargs()
  print "trying to remove functions %s from %s" % (funcls, fninp)
  s = open(fninp).readlines()
  s = rmfuncs(s, funcls, force)
  if fnout:
    print "writing", fnout
    open(fnout, "w").writelines(s)
  else:
    print ''.join(s)


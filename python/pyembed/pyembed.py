#!/usr/bin/env python

''' embed python '''

import os, sys, re, getopt



verbose = 0
doall = False



def getfuncs(s):
  ''' collect functions '''

  ls = []
  for i in range(len(s)):
    ln = s[i].strip()
    m = re.match(r"\s*def\s+(\w+?)\s*?\(.*", ln)
    if m: ls += [ m.group(1), ]
  return ls



def getmodls(s):
  ''' collect a list of imported modules '''

  ls = []
  for i in range(len(s)):
    ln = s[i].strip()
    m = re.match(r"\s*import\s+(.*?)\s*;?\s*$", ln)
    if not m:
      m = re.match(r"\s*from\s+(.*?)\s+import\s+.*$", ln)
    if not m: continue
    for a in m.group(1).strip().split(","):
      mod = a.strip()
      # exclude commonly used modules, which need no embedding
      if mod in ("os", "sys", "re", "string", "getopt", "shutil",
          "subprocess", "glob", "Gnuplot",
          "struct", "difflib", "datetime", "calendar", "heapq",
          "collections", "bisect", "array", "sets", "mutex",
          "numbers", "cmath", "math", "fractions", "random",
          "pickle", "zlib", "gzip", "bz2", "zipfile", "tarfile"):
        continue
      ls += [ mod, ]
  return ls



def embed0(src, mod, imported):
  """ embed `mod' into src """

  # prepare the module to be inserted
  fnmod = mod + ".py"
  if not os.path.exists(fnmod):
    print "cannot import %s" % fnmod
    return src
  srcmod = open(fnmod).readlines()
  if len(srcmod) and srcmod[0].startswith("#!"):
    srcmod = srcmod[1:]
  srcmod  = ["\n\n################### embedded module '%s' begins ####################\n" % mod, ] + srcmod
  srcmod += ["################### embedded module '%s' ends ######################\n\n\n" % mod, ]

  first = -1
  # 1. remove all statements of importing `mod'
  for i in range(len(src)):
    ln = src[i].rstrip()

    # remove comments
    k = ln.find("#")
    if k >= 0: ln = ln[:k]

    # import patterns
    repls = [
      [r"(\s*?)import\s+%s\s*;?\s*$" % mod, ""],
      [r"(\s*?)import\s+%s\s*,(.*)$" % mod, r"\1import \2" + "\n"],
      [r"(\s*?)import\s+(.*),\s*%s\s*;?\s*$" % mod, r"\1import \2" + "\n"],
      [r"(\s*?)import\s+(.*),\s*%s\s*,(.*)$" % mod, r"\1import \2, \3" + "\n"],
      [r"(\s*?)from\s+%s\s+import\s+.*$$" % mod, ""],
    ]

    # matching patterns
    for (find, repl) in repls:
      m = re.match(find, ln)
      if m:
        src[i] = re.sub(find, repl, ln)
        #print "find %s, repl %s\nline %s: %s" % (find, repl, i, ln)
        # the insertion point must have no indent
        if first < 0 and not m.group(1): first = i
        break
    else:
      continue

  # 2. insert the source module
  if not mod in imported:
    # try to find a suitable start point of insertion
    if first < 0:
      first = 0
      while first < len(src): # use the first blank line
        if not src[first].strip(): break
        first += 1
      while first < len(src): # now skip over a blank block
        if src[first].strip(): break
        first += 1
      # also skip a description block
      bstr = None
      if src[first].startswith("'''"): bstr = "'''"
      if src[first].startswith('"""'): bstr = '"""'
      if bstr:
        if bstr in src[first][3:]:
          # a single line description: ''' blah blah blah... '''
          first += 1
        else:
          # a multiple line description
          first += 1
          while first < len(src):
            if bstr in src[first]: break
            first += 1
          else: # retreat one line, for we have reach the end
            first -= 1
          first += 1

        # now skip over a blank block
        while first < len(src):
          if src[first].strip(): break
          first += 1

    if verbose:
      print "module '%s' is to be inserted before line %s" % (mod, first)
    src = src[:first] + srcmod + src[first:]

  # 3. change function references
  for i in range(len(src)):
    src[i] = re.sub(r"(?<!\w)%s\." % mod, "", src[i])

  imported[ mod ] = True
  return src



def embed(src, modls):
  ''' embed a list of modules '''
  imported = {} # dictionary of imported modules
  # iterate the module-insertion process, because modules may import each other
  for i in range(len(modls) + 1):
    if verbose: print "starting iteration %d" % i
    src0 = src[:]
    for mod in modls:
      src = embed0(src, mod, imported)
    if ''.join(src0) == ''.join(src):
      if verbose: print "finishing at iteration %d" % i
      break
  return src



def usage():
  print """  Usage:
    %s [OPTIONS] input.py
  Embed the contents of some modules into `input.py'
  Copyright (c) 2013 Cheng Zhang

  OPTIONS:
  -m, --mod=        list all modules to be embedded, separated by `,'
  -a, --all         embed all modules
  -o, --output=     output file
  -v                increase the verbosity
  --verbose=        set the verbosity
  -h, --help        help
  """ % sys.argv[0]
  exit(1)



def doargs():
  ''' handle the command-line arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "am:o:vh",
         ["all", "mod=", "output=", "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err)
    usage()

  global doall, verbose

  verbose = 0
  fninp = fnout = None
  modls = []

  for o, a in opts:
    if o in ("-a", "--all",):
      doall = True
    elif o in ("-m", "--mod",):
      modls = [ s.strip() for s in a.strip().split(",") ]
    elif o in ("-o", "--output",):
      fnout = a
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

  if not modls and doall:
    modls = getmodls(open(fninp).readlines())

  return fninp, fnout, modls



if __name__ == "__main__":
  fninp, fnout, modls = doargs()
  print "embedding %s into %s" % (modls, fninp)
  s = open(fninp).readlines()
  s = embed(s, modls)
  if fnout:
    print "writing", fnout
    open(fnout, "w").writelines(s)
  else:
    print ''.join(s)


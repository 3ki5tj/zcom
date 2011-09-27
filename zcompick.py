#!/usr/bin/env python

''' create a specialized version from zcom.h '''

import os, sys, re, shutil, getopt, filecmp

verbose = 0
fninput = "zcom.h"
fnoutput = "zcom1.h"
keys = ['def', 'util', 'rng', 'cfg', 'trace', 'endn', 'ss']
prefix = "ZCOM_"

def mksmall(input, output, goodkeys):
  ''' created a specialized version 'output' (zcom1.h) from 'input' (zcom.h) 
      with a few modules, given by 'goodkeys' '''

  oldsrc = open(input, 'r').readlines()
  goodkeys = transkeys( goodkeys, prefix )
  badkeys = getbadkeys(oldsrc, goodkeys)
  newsrc = []
  
  lnum = 1
  level = 0    # the #if/#ifdef/#ifndef block level
  ignlev = -1  # if >= 0, the current line is ignored
  toclev = -1

  def haskey(ln, keys):
    ''' if a line ln has any of keys '''
    ln = ln.strip()
    for key in keys:
      if re.search("\s" + key + "$", ln):
        if verbose >= 1: print "IGNOR + %s, #%s: %s" % (level, lnum, lin),
        return 1
    return 0

  for line in oldsrc: # do line by line
    lin = line.lstrip()

    # increase the level
    if lin.startswith("#if"):
      level += 1
      if verbose >= 2: print "LEVEL + %s, #%s: %s" % (level, lnum, lin),

      # the "table of contents" led by ZCOM_PICK
      if toclev < 0:
        if lin.startswith("#ifndef ZCOM_PICK"):
          toclev = level
      elif ignlev < 0 and lin.startswith("#ifndef") and haskey(lin, badkeys):
        ignlev = level

      # ignore
      if ignlev < 0:
        if level == 1 and lin.startswith("#ifdef"):
          if haskey(lin, badkeys): ignlev = level
        elif lin.find("_LEGACY") >= 0:
          ignlev = level

    if ignlev < 0: newsrc += [line]

    # decrease the level
    if lin.startswith("#endif"):
      if toclev >= 0 and toclev == level:
        toclev = -1

      if ignlev >= 0 and ignlev == level:
        if verbose >= 1: print "IGNOR - %s, #%s: %s" % (level, lnum, lin),
        ignlev = -1

      level -= 1
      if verbose >= 2: print "LEVEL - %s, #%s: %s" % (level, lnum, lin),
      if level < 0:
        print "negative level %s, #%s: %s" % (level, lnum, lin)
        raise Exception

    lnum += 1 # increase the line number
  
  newsrc = rmblank(newsrc, 2) # remove more than 2 successive blank lines 
  writefile(output, newsrc, goodkeys, badkeys)

def getbadkeys(lines, excls):
  ''' search all keys in the toc section, and exclude good keys `excl' 
      return a list of bad keys '''
  toclev = -1
  level = 0
  list = []
  lnum = 1
  for line in lines:
    lin = line.strip()
    if lin.startswith("#if"):
      level += 1

    if lin.startswith("#ifndef"):
      if level == 1 and (toclev < 0) and lin.find("ZCOM_PICK") >= 0:
        toclev = level
      elif toclev >= 0 and level == 2:
        key = lin.split()[1]
        if not (key in excls):
          list += [ key ]

    if lin.startswith("#endif"):
      if toclev >= 0 and toclev == level:
        toclev = -1
        return list
      level -= 1
    lnum += 1
  return list

def rmblank(src, tol):
  ''' remove more than tol successive blank lines '''
  i = 0
  last = -1 # last nonblank line
  newsrc = []
  for line in src:
    lin = line.strip()
    if lin != "":
      last = i
    if i - last <= 2:
      newsrc += [line]
    i += 1
  return newsrc

def writefile(fn, src, goodkeys, badkeys):
  ''' src to file fn '''
  # write to a temporary file first then sync
  fntmp = fn + ".tmp"
  open(fntmp, 'w').write(''.join(src))
  if not os.path.exists(fn) or not filecmp.cmp(fntmp, fn):
    print "good: %s\n bad : %s" % (goodkeys, badkeys)
    shutil.copy2(fntmp, fn)
    print "update %s" % fn
    ret = 1
  else:
    print "no need to update %s" % fn
    ret = 0
  os.remove(fntmp)
  return ret

def transkeys(keys, pfx):
  ''' add prefix to keys '''
  newkeys = []
  for k in keys:
    key = k.upper()
    if not key.startswith(pfx):
      key = pfx + key
    newkeys += [ key ]
  return newkeys

def usage():
  ''' print usage and die '''
  print sys.argv[0], "[Options]"
  print " -i, --input=\n\tfollowed by input file, default %s" % fninput
  print " -o, --output=\n\tfollowed by output file, default %s" % fnoutput
  print ' -k, --keys=\n\tfollowed by a list of keys, default "%s"' % ' '.join(keys)
  print " -p, --prefix=\n\tfollowed by the prefix, default %s" % prefix
  print " -v, --verbose=\n\tfollowed by verbose level"
  sys.exit(1)

def doargs():
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hi:o:k:v:",
         ["help", "input=", "output=", "keys=", "key=", "verbose="])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global fninput, fnoutput, keys, prefix, verbose

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-i", "--input"):
      fninput = a
    elif o in ("-o", "--output"):
      fnoutput = a
    elif o in ("-k", "--key", "--keys"):
      keys = a.split()
    elif o in ("-p", "--prefix"):
      prefix = a
    elif o in ("-h", "--help"):
      usage()
    else:
      print "unknown option '%s'" % o
      usage()

if __name__ == "__main__":
  doargs()
  mksmall(fninput, fnoutput, keys)


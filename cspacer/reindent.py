#!/usr/bin/env python

''' format the leading spaces and tabs
    Copyright (c) 2013 Cheng Zhang

    Main functions:
    * guessindent():    guess the indent size
    * tab2sp():         change leading tabs to spaces
    * reindent():       main function, reformat the leading indents
      reindentf():      wrapper for reindent, accepts a file

    Helper functions:
    * cleanline():      remove comments, then strip()
    * getindent():      get the leading indent of the line
    * switchpp():       switch the current preprocessor state
    * switchcmt():      switch the current comment state
    * gcd():            greatest common divisors
'''

import os, sys, getopt, re
import fileglob, cfmt

verbose = 0
defindent = "  "
forcedef = False
backup = True
fnout = None
emacstag = True  # trust Emacs tag for tab size

strict = False
maxind = 80 # maximial number of leading characters

# flags for preprocessor and comment states
PPFLAG = 1
CMTFLAG = 2



''' ################## Helper functions begins ##################### '''

def cleanline(ln):
  ''' strip away comments and spaces '''

  ln = re.sub(r"//.*", "", ln).strip()
  ln = re.sub(r"/\*.*?\*/", "", ln).strip()
  # filter out a partial comment `/* .... '
  ln = re.sub(r"/\*.*?$", "", ln).strip()
  ln = re.sub(r"^.*?\*/", "", ln).strip()
  return ln



def getindent(ln):
  ''' get the leading indent '''

  for k in range(len(ln)):
    if not ln[k].isspace(): return ln[:k]
  else:
    return ln



def switchpp(ln, pplevel, inpp, incmt):
  ''' switch the current state for preprocessors '''

  if not incmt:
    if ln.startswith("#"):
      inpp = True
      if ln.startswith("#if"):
        pplevel += 1
      elif ln.startswith("#endif"):
        pplevel = max(pplevel - 1, 0)

  inppthis = inpp

  # compute the pp state of the next line
  if inpp:
    # end the processor, unless there's a line continuation
    if not ln.endswith("\\"):
      inpp = False

  return pplevel, inpp, inppthis



"""
def getprevln(s0, i, inppcmt = None):
  ''' get the index of the preceeding code line before line `i' '''

  iprev = max(i - 1, 0)
  while iprev > 0:
    pln = s0[iprev].strip()
    pgood = True

    # always skip a blank line
    if pln == "": pgood = False
    
    # the previous line is a label, e.g., EXIT:
    if re.match("\w+:", pln): pgood = False
    
    # comment or preprocessor 
    if inppcmt and inppcmt[iprev]: pgood = False
    
    if pgood: break
    iprev -= 1
  return iprev
"""


def gcd(a, b):
  ''' greatest common divisor '''
  if a > b: a,b = b,a
  if a == 0: return b
  else: return b % a


''' ################## Helper functions ends ##################### '''






def guessindent(lines):
  ''' guess the basic unit of leading tabs '''
  unitind = defindent

  # 1. if the user insist
  if forcedef: return defindent

  # 2. try to read the `tab-width' tag
  if emacstag:
    for i in range(len(lines)):
      if len(lines[i].strip()): break
    else:
      return defindent

    ln = lines[i]
    if ln.startswith("/*"):
      m = re.search("tab-width:\s*([0-9]+);", ln)
      if m:
        if verbose:
          print "Tab-size according to `tab-width':", m.group(1)
        ntab = int(m.group(1))
        return " " * ntab

  # 3. try to guess the indent size
  lntab = 0 # number of lines that uses tabs
  lnsp = [0] * maxind # number of lines that uses spaces
  pplevel = 0
  inpp = False  # in preprocessors
  incmt = False # in comments
  lnum = 0
  for line in lines:
    lnum += 1
    ln = line.rstrip()

    # preprocessor control
    pplevel, inpp, inppthis = cfmt.switchpp(ln, pplevel, inpp, incmt)
    if inppthis: continue

    # comment control
    incmt, incmtthis = cfmt.switchcmt(ln, incmt)
    if incmtthis: continue
    
    ind = getindent(ln)
    if len(ind) == 0: continue
    
    # update the # of lines that use Tab as indents
    l = ind.find("\t")
    if l >= 0:
      ind = ind[:l]
      lntab += 1
    
    # compute the number of spaces in the indent
    ns = len(ind)
    if ns == 0: continue
    if ns % 2 != 0 and verbose >= 2:
      print "strange indents", ns, ":", lnum, ":", line.strip()
    if ns < len(lnsp): lnsp[ns] += 1


  nsp = lnsp.index( max(lnsp) )
  if nsp % 2 == 0: # use a heuristic formula
    nsp = 8
    if lnsp[4] >= 0.2 * lnsp[8]:
      nsp = 4
      if lnsp[2] >= 0.2 * lnsp[4]:
        nsp = 2
  if verbose:
    print "tabs:", lntab, "lines, spaces:", lnsp[1:9], "lines, indent", nsp

  if lntab > lnsp: # more lines uses tabs
    return "\t"
  else: # more lines uses spaces
    return " " * nsp



def tab2sp(s0):
  ''' remove leading tabs in the indent to spaces
      called in `reindent'  '''

  # guess the soft tab size
  unitind = guessindent(s0)

  s1 = [ln.rstrip() + '\n' for ln in s0]
  n = len(s0)
  i = 0
  while i < n:
    ln = s1[i].rstrip()
    lnstrip = ln.strip()
    i += 1

  return s0[:]



def reindent(s0):
  ''' reformat the leading indents
      when `{' is encountered, an indent is added on the next line
      when `}' is encountered, the indent is removed
      for other delibrate indents, 
      the algorithm looks for the incremental indent difference
      between the current line and previous line, and try to
      the mimic the difference in the output '''

  # remove trailing spaces
  s0 = [ln.rstrip() + '\n' for ln in s0]
  # change leading tabs to spaces
  s1 = tab2sp(s0)

  unitind = guessindent(s0)
  level = 0
  levels = [ level ] * len(s0)
  indents = [ "" ] * 100

  iprev = 0
  pplevel = 0
  inpp = False  # in preprocessor
  incmt = False # in comments
  inppcmt = [0] * len(s0)
  n = len(s0)

  for i in range(n):
    ln = s0[i].rstrip()
    lnstrip = ln.strip()
    if len(lnstrip) == 0: # blank line
      continue

    # preprocessor control
    pplevel, inpp, inppthis = cfmt.switchpp(lnstrip, pplevel, inpp, incmt)
    if inppthis:
      inppcmt[i] = PPFLAG
      continue

    # control comments
    incmt, incmtthis = cfmt.switchcmt(lnstrip, incmt)
    if incmtthis:
      inppcmt[i] += CMTFLAG
      if i > 0: levels[i] = levels[iprev]
      continue

    # start formating
    lnclean = cleanline(ln)

    follow = 0
    indincr = inddecr = ""
    # if this line has no indent, reset
    if not ln[0].isspace():
      # avoid resetting if it is a tag line like `EXIT:'
      if re.match(r"\w+:", lnclean):
        s1[i] = s0[i].strip() + '\n'
        continue
      else:
        level = 0
    elif i > 0:
      # compute the indent difference from the previous line
      indthis = getindent(s0[i])
      indprev = getindent(s0[iprev])
      # follow this line if it has the same indent as the previous one
      if indthis == indprev:
        follow = 1
      # register the amount of indent in the source
      if indthis.startswith(indprev):
        indincr = indthis[len(indprev):]
      elif indprev.startswith(indthis):
        inddecr = indprev[len(indthis):]

    # adjust the indent level, before the current line
    if lnclean.startswith("}"):
      level = max(0, level - 1)
    levels[i] = level

    # compute the indent for the current line 
    if follow: # following the previous line's indent
      indent = getindent(s1[iprev])
    else:
      if strict:
        indent = unitind * level
      else:
        indent = indents[level]
      # if this line shares the same level with the previous one
      # but it appears to have longer/shorter indent,
      # then try to mimic the source
      if i > 0 and level == levels[iprev]:
        indprev = getindent( s1[iprev] )
        if len(indincr) > 0:
          indent = indprev + indincr
        elif len(inddecr) > 0 and indprev.endswith(inddecr):
          indent = indprev[:-len(inddecr)]

    # apply the indent for the current line
    s1[i] = indent + s0[i].strip() + '\n'

    # DEBUGGING BLOCK
    if verbose == 10:
      print "iprev %s, i %s, [%s,%s]" % (iprev, i, getindent(s0[iprev]), getindent(s0[i]))
      print "indent [%s], level %d, %d, indprev [%s] indincr [%s] inddecr [%s] follow %s, inppcmt %s\n%s" % (
          indent, level, levels[iprev], indprev, indincr, inddecr, follow, inppcmt[i], lnstrip)
      raw_input()

    # adjust the indent level for the following lines
    iref = -1
    if lnclean.endswith("{"):
      iref = i
      # if we encounter an ending `) {'
      # try to find its starter if/else/for/while
      if re.search(r"\)\s*{$", lnclean):
        while iref >= 0:
          lnref = cleanline( s0[iref] )
          if lnref == "":  # gone too far
            iref += 1
            break
          if ( lnref.startswith("if")
            or lnref.startswith("for")
            or lnref.startswith("while")
            or re.search(r"^(})?\s*else", lnref) ):
            #print "ref-ln %d -> %d" % (i, iref)
            break
          iref -= 1
      ss = getindent(s1[iref])
      # try to update also the indent of this level
      # as long as this is not the ground level,
      if level > 0: indents[ level ] = ss
      indents[ level+1 ] = ss + unitind
      level += 1

    if s1[i] != s0[i] and verbose >= 2:
      print "%s, prev-level %s, level %s, next-level %s, indent %s, [%s], incmt %s, iprev %s, iref %s\nold:%snew:%s" % (
          i, levels[iprev], levels[i], level, len(indent), indent, incmt, iprev, iref,
          s0[i], s1[i]),
      if verbose >= 3:
        raw_input("indent %s, follow %s, [%s], [%s]" % (
          indent, follow, getindent(s0[iprev]), indents[:level+2] ) )
    
    # register the last good code line
    iprev = i
  # end of the main loop
  
  return s1



def reindentf(fn):
  ''' reindent a file, wrapper '''

  try:
    s0 = open(fn).readlines()
  except Exception:
    print "cannot open %s" % fn
    return

  s1 = reindent(s0)
 
  # be careful when rewritting the file
  if ''.join(s0) == ''.join(s1):
    print "keeping", fn
  else:
    global fnout
    if not fnout:
      fnout = fn
      if verbose >= 1:
        yesno = raw_input("ok?").strip()
        if not yesno or not yesno[0] in "yYoO":
          print "abort"
          exit(1)
      if backup:
        for i in range(1, 100):
          fnbak = fn + ".bak" + str(i)
          if not os.path.exists(fnbak): break
        open(fnbak, "w").writelines(s0)
    print "writing", fnout
    open(fnout, "w").writelines(s1)



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] input"
  print """
  Reformat leading indents (spaces and tabs)

  OPTIONS:

   -t:                by default, use tab instead of spaces
   -s n, --spaces=n:  the default number `n' of spaces of an indent
   -f, --force:       always use the default indent
   -l, --links:       include symbolic links
   -R, --recursive:   apply to subdirectories, if `input' is
                      a wildcard pattern like `*.c',
                      the pattern must be quoted as '*.c'
   -o, --output:      output file
   --noemacs:         don't trust Emacs tag for the tab size
   -v:                be verbose
   --verbose=n:       specify the verbosity
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "ts:flRo:hv",
         ["spaces=", "force", "links", "recursive", "output=",
          "noemacs", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, defindent, forcedef, fnout, emacstag

  recur = links = False

  for o, a in opts:
    if o in ("-f", "--force",):
      forcedef = True
    elif o in ("-s", "--spaces",):
      defindent = " " * int(a)
    elif o in ("-t", "--tab",):
      defindent = "\t"
    elif o in ("-R", "--recursive",):
      recur = True
    elif o in ("-l", "--links",):
      links = True
    elif o in ("-o", "--output",):
      fnout = a
    elif o in ("--noemacs",):
      emacstag = False
    elif o in ("-v",):
      verbose = 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  # common C code
  pats = '''*.c *.cpp *.h *.hpp *.java'''.split()
  if len(args) > 0:
    # parse the pattern in each argument
    pats = [ a for pat in args for a in pat.split() ]

  # compile a list of files
  ls = fileglob.fileglob(pats, links, recur)
  if len(ls) <= 0: print "no file for %s" % pats
  return ls



def main():
  fns = doargs()
  for fn in fns: reindentf(fn)


if __name__ == "__main__":
  main()


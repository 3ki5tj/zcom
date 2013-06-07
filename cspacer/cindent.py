#!/usr/bin/env python

''' format the leading spaces and tabs
    Copyright (c) 2013 Cheng Zhang

    Main functions:
    * tab2sp():         change leading tabs to spaces
    * reindent():       main function, reformat the leading indents
      reindentf():      wrapper for reindent, accepts a file
    * guessindent():    guess the indent size

    Helper class CParser
      * switchpp():       switch the current preprocessor state
      * switchcmt():      switch the current comment state

    Helper functions:
    * cleanline():      remove comments, then strip()
    * getindent():      get the leading indent of the line
    * gcd():            greatest common divisors
    * endingcmt():      if a code line starts a block comment
'''

import os, sys, getopt, re

verbose = 0
defindent = "  "
forcedef = False
backup = True
fnout = None
emacstag = True  # trust Emacs tag for tab size
fmtcmt = True # treat comment lines as code

strict = False
maxind = 80 # maximial number of leading characters




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



def gcd(a, b):
  ''' greatest common divisor '''
  if a > b: a,b = b,a
  if a == 0: return b
  else: return b % a



def endingcmt(ln, incmt = False):
  ''' check if this line ends with a block comment
        a=b; /* hahaha */ b=3; /* comment
      NOTE: it will fail if `/*' or '*/' is hidden in a string '''

  if incmt: # finish the current comment
    m = re.match(r".*?\*/(.*)", ln)
    if m: ln = m.group(1)
    else: return True # the block comment remains

  while 1:
    m = re.match(r".*?(/\*.*?\*/)(.*)", ln)
    m1 = re.match(r".*?(/\*).*", ln)
    m2 = re.match(r".*?(//).*", ln)

    if m2: # has `//'
      if not m1: return False
      # both '//' and '/*' exist
      if m2.start(1) < m1.start(1):
        # `//' precedes `/*'
        #print "line comment:", ln[m2.start(1):]
        return False
      elif not m:
        # `/*' precedes `//', but no '*/'
        return True
      else: pass # closed `/* ... */', pass through

    if not m: # no more completed `/* ... */'
      return (m1 != None)
    else:
      ln = m.group(2)

''' ################## Helper functions ends ##################### '''



''' ############## Helper class CParser begins ################### '''
class CParser:
  ppn = 0  # preprocessor level
  inpp = False  # currently in a preprocessor
  incmt = False

  def __init__(self, ppn = 0, inpp = False, incmt = False):
    self.ppn = ppn
    self.inpp = inpp
    self.incmt = incmt


  def switchpp(self, ln):
    ''' switch the current state for preprocessors
        according to the input line `ln'
        return if the current line is in a preprocessor block '''

    ln = ln.strip()
    # ignore the preprocessor, if it is in a comment
    if not self.incmt:
      if ln.startswith("#"):
        self.inpp = True
        if ln.startswith("#if"):
          self.ppn += 1
        elif ln.startswith("#endif"):
          self.ppn = max(self.ppn - 1, 0)

    inppthis = self.inpp

    # compute the pp state of the next line
    if self.inpp:
      # end the processor, unless there's a line continuation
      if not ln.endswith("\\"):
        self.inpp = False
    return inppthis


  def switchcmt(self, ln):
    ''' switch the current state for block comments
        according to the input line `ln'
        an embedded comment does not counts:
          a = 3; /* assignment */ b = 4;
        return if the current line is in a comment block '''

    ln = ln.strip()
    incmtthis = self.incmt
    if not self.incmt and ln.startswith("/*"):
      incmtthis = True
    self.incmt = endingcmt(ln, self.incmt)
    return incmtthis

"""
def printppcmt(s):
  ''' test: print preprocessor or comment lines '''
  cp = CParser()
  for i in range(len(s)):
    inpp = cp.switchpp(s[i])
    incmt = cp.switchcmt(s[i])
    if inpp or incmt: print "%5d:%d,%5s|%5s|%s" % (i+1, cp.ppn, cp.inpp, cp.incmt, s[i]),
"""

''' ############## Helper class CParser ends ################### '''




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
  cp = CParser()
  lnum = 0
  for line in lines:
    lnum += 1
    ln = line.rstrip()

    # skip preprocessor and comments
    if cp.switchpp(ln): continue
    if cp.switchcmt(ln): continue

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
    if ns % 2 != 0 and verbose >= 3:
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



def lntab2sp(s, tab):
  ''' convert tab into spaces '''
  s1 = ""
  for k in range(len(s)):
    if s[k] == '\t':
      s1 += " " * int( (len(s1) + tab) / tab ) * tab
  return s1



def tab2sp(s0, tabsize = 4):
  ''' remove leading tabs in the indent to spaces
      called in `reindent'  '''

  s1 = [ln.rstrip() + '\n' for ln in s0]
  n = len(s0)
  iprev = 0
  cp = CParser()
  for i in range(n):
    # skip preprocessor and comments
    if cp.switchpp(s0[i]): continue
    if cp.switchcmt(s0[i]): continue

    ind0 = getindent(s0[i])

    # skip if the previous line has no tab
    if '\t' in ind0:
      # compute the proper tabsize
      indprev0 = getindent(s0[iprev])
      indprev1 = getindent(s1[iprev])
      if indprev0 == ind0:
        # if the indent is the same as the previous line, use the old tab
        ind1 = indprev1
      else:
        # try different tab sizes, and choose the proper one
        for tb in [2, 4, 8]:
          ind1 = lntab2sp(ind0, tb)
          if len(ind1) >= len(indprev1):
            break
        if verbose >= 2:
          print "TAB2SP line %s, [%s|%s]%s" % (i+1, indprev1, ind1, s0[i].strip())
      # apply the new indent
      s1[i] = ind1 + s1[i].lstrip()
    iprev = i
  return s1



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
  unitind = guessindent(s0)
  # change leading tabs to spaces
  if unitind != '\t':
    s1 = tab2sp(s0, len(unitind))
  else:
    s1 = s0[:]

  level = 0
  levels = [ level ] * len(s0)
  indents = [ "" ] * 100

  iprev = 0
  cp = CParser()
  n = len(s0)

  for i in range(n):
    ln = s0[i].rstrip()
    lnstrip = ln.strip()
    if len(lnstrip) == 0: # blank line
      continue

    # preprocessor control
    if cp.switchpp(lnstrip): continue

    # control comments
    incmt = cp.switchcmt(lnstrip)
    if not fmtcmt: # skip comments
      if incmt: continue

    # start formating
    lnclean = cleanline(ln)

    follow = 0
    indincr = inddecr = ""
    # if this line has no indent, reset
    if not ln[0].isspace():
      # avoid resetting if it is a label line like `EXIT:'
      if re.match(r"\w+:", lnclean):
        s1[i] = s0[i].strip() + '\n'
        continue
      elif level != 0:
        if incmt: # a quick fix for mdrun.c
          continue
        if verbose >= 3:
          print "RESET indent level %s at line %s: %s" % (level, i+1, lnstrip)
        # reset the level, nasty trick, should be avoided
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
      if verbose >= 3:
        print "UNDENT line", i+1, "level", level, ":", s1[i].strip()
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
    if 0 and i == 381:
      print "DEBUG iprev %s, i %s, [%s,%s]" % (iprev, i, getindent(s0[iprev]), getindent(s0[i]))
      print "DEBUG indent [%s], level %d, %d, indprev [%s] indincr [%s] inddecr [%s] follow %s\n%s" % (
          indent, level, levels[iprev], indprev, indincr, inddecr, follow, lnstrip)
      raw_input()

    # adjust the indent level for the ensuing lines
    iref = -1
    if lnclean.endswith("{"):
      iref = i
      # if we encounter an ending `) {'
      # try to find its starter if/else/for/while
      # s.t. we can reset the indent level
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
      if verbose >= 3:
        print "INDENT line", i+1, "level", level, "iref", iref+1, ":", s1[iref].strip()
      # try to update also the indent of this level
      # as long as this is not the ground level,
      if level > 0: indents[ level ] = ss
      indents[ level+1 ] = ss + unitind
      level += 1

    # The difference does not count the trailing spaces
    if s1[i] != s0[i] and verbose >= 2:
      print "DIFF: %s, prev-level %s, level %s, next-level %s, indent %s, [%s], incmt %s, iprev %s, iref %s\nold:%snew:%s" % (
          i+1, levels[iprev], levels[i], level, len(indent), indent, cp.incmt, iprev, iref,
          s0[i], s1[i]),

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
   -R, --recursive:   apply to subdirectories, if `input' is
                      a wildcard pattern like `*.c',
                      the pattern must be quoted as '*.c'
   -L, --nolinks:     skip symbolic links
   -o, --output:      output file
   --noemacs:         don't trust Emacs tag for the tab size
   -v:                be verbose
   --verbose=n:       specify the verbosity
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "ts:fRLo:hv",
         [ "spaces=", "force", "recursive",
           "nolinks", "output=",
           "noemacs", "verbose=", "help", ] )
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, defindent, forcedef, fnout, emacstag

  recur = False
  links = True

  for o, a in opts:
    if o in ("-f", "--force",):
      forcedef = True
    elif o in ("-s", "--spaces",):
      defindent = " " * int(a)
    elif o in ("-t", "--tab",):
      defindent = "\t"
    elif o in ("-R", "--recursive",):
      recur = True
    elif o in ("-L", "--nolinks",):
      links = False
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

  ls = args
  try: # limit the dependence on argsglob
    from zcom import argsglob
    ls = argsglob(args, "*.c *.cpp *.h *.hpp *.java", recur = recur, links = links)
  except ImportError: pass
  if len(ls) <= 0: print "no file for %s" % args
  return ls



def main():
  fns = doargs()
  for fn in fns: reindentf(fn)


if __name__ == "__main__":
  main()


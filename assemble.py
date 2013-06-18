#!/usr/bin/env python

'''
integrate target.c and target.h into fn_host
debug code in target.c is removed first by calling rmdbg
C++ compments are also removed
'''

import os, sys, inspect, shutil, getopt, re

# import functions from the subdirectory `python'
from python import rmdbg, depls


# module attributes, set before use
fnhost      = "zcom.h"
fnoutp      = None
strcls      = "STRCLS"
prefix      = "ZCOM_"
verbose     = 1


def strip_def(src, fn):
  ''' strip away the
    #ifndef FILE__
    #define FILE__
    #endif
  triplet, and everything outside '''

  # form an upper case `tag' from the file name
  tag = fn.replace(".", "_").upper() + "__"

  n = len(src)

  # 1. find `#ifndef'
  for i in range(n):
    ln = src[i].strip()
    if ln.startswith("#ifndef") and ln[7:].strip().startswith(tag):
      if verbose > 1: print "#ifndef of %s found in line %s: %s" % (fn, i, ln)
      break
  else:
    print "no #ifndef in %s" % fn
    raise Exception
    return src # no "#ifdef" found

  # 2. look for the ensuing `#define'
  start = i + 1
  if not src[start].lstrip().startswith("#define"):
    return src
  if verbose > 1:
    print "#define found in line", start

  start += 1  # beginning of the main code
  # 3. search for the final "#endif"
  for i in range(n - 1, start, -1):
    if src[i].lstrip().startswith("#endif"):
      end = i
      break
  else:
    print "no #endif in %s" % fn
    raise Exception

  return src[start : end]



def get_includes(src):
  ''' get a list of files included in src '''

  src = [ln for ln in src if ln.strip().startswith("#include")]
  fns = []
  for ln in src:
    m = re.search('\s*#include\s*"(.*?)"', ln)
    if m: fns += [ m.group(1) ]
  return fns



def add_storage_class(hdr, prefix):
  ''' add storage class `prefix' to global functions and variables in the header '''

  for i in range(len(hdr)):
    line = hdr[i]
    arr = line.split()
    if len(arr) > 0:
      if line[0:1].isspace():
        continue
      else:
        first_word = arr[0]
    else:
      continue

    #  NOTE: we cannot handle more complex cases
    if first_word in ("void", "char", "int", "unsigned", "long",
        "real", "float", "double", "const",) or first_word.endswith("_t"):
      hdr[i] = prefix + " " + hdr[i]
  return hdr



def remove_cpp_comments(src):
  ''' remove C++ style comments
      only a line that starts with // is included '''

  i = 0
  while i < len(src):
    s = src[i].strip()
    if s.startswith("//"):
      if verbose >= 2:
        print "Remove comment [%s]" % s; raw_input()
      src = src[:i]+src[i+1:]
    else:
      i += 1
  return src



def insert_module(src, name, smodule):
  ''' find the block specified by
        #ifdef  name
        #ifndef name__
        #define name__
        ...
        #endif
        #endif
      in `src' and insert `smodule' into `...'  '''

  n = len(src)
  namein = name + "__"
  plevel = 0
  bstart = -1
  bend = -1

  for i in range(n):
    if (plevel == 0 and bstart < 0 and
        src[i].startswith("#ifdef") and
        src[i][6:].strip().startswith(name) and
        src[i+1].startswith("#ifndef") and
        src[i+1][7:].strip().startswith(namein) and
        src[i+2].startswith("#define") and
        src[i+2][7:].strip().startswith(namein) ):
      bstart = i + 3
      if verbose > 1:
        print "out loop starts from:", i, "line:", src[i].rstrip()

    if (plevel == 2 and
        bstart > 0 and
        src[i].startswith("#endif") and
        src[i+1].startswith("#endif") ):
      bend = i
      if verbose > 1:
        print "out loop ends at:", i+1, "line:", src[i+1].rstrip()
      break

    # update #if/#endif loops
    if src[i].strip().startswith("#if"):
      plevel += 1
    elif src[i].strip().startswith("#endif"):
      plevel -= 1

  if bstart > 0 and bend > 0:
    return src[:bstart] + smodule + src[bend:]
  else:
    print "cannot find the desired module", name
    return src



def builddeps(srcls):
  ''' build dependencies '''
  n = len(srcls)
  #for i in range(n): print "%d: %s" % (i, srcls[i][0])

  # 1. build dependencies
  deps = [""]*n
  depsX = [""]*n
  for i in range(n):
    mod = srcls[i][0]
    deps[i] = []
    depsX[i] = []
    for j in range(n):
      if j == i: continue
      d = srcls[j][0]
      file_h = os.path.join(mod, d+'.h')
      if os.path.exists(file_h):
        deps[i] += [j,]
      file1_h = os.path.join(mod, "_"+d+'.h') # virtual dependence
      file1_c = os.path.join(mod, "_"+d+'.c') # virtual dependence
      if os.path.exists(file_h) or os.path.exists(file1_h) or os.path.exists(file1_c):
        depsX[i] += [j,]
    #if (len(deps[i])): print mod, deps[i]

  # 2. sort modules according to dependencies
  # d[0], d[1], ..., d[n-1] is a permutation of n
  # d[0] is most basic, i.e., no dependencies
  # d[k] can depend on (d[0], d[1], ... d[k-1])
  d, cycl = depls.depsort(depsX)
  if cycl:
    print "cyclic dependency detected"
    for i in cycl:
      print "%d (%s) -> " % (i, srcls[i][0]),
    raise Exception

  # 3. build a new sorted src list and dependence list
  newsrcls = [0]*n
  newdeps = [0]*n
  for i in range(n):
    id = d[i]
    newsrcls[i] = srcls[id]
    newdeps[i] = sorted([d.index(j) for j in deps[id]])

  # 4. prune and minimize the dependence list
  # if C depends on A, B, B depends on A
  # it can be simplified as
  # C on B, B on A
  mindeps = [0]*n # a minimal dependent
  for i in range(n):
    m = len(newdeps[i])
    depi = []
    for j in range(m):
      needj = 0
      dj = newdeps[i][j]
      for k in range(j+1, m):
        dk = newdeps[i][k]
        # if k depends on j too, then no need to put j in the list
        if dj in newdeps[dk]:
          break
      else: # j is necessary
        needj = 1
      if needj: depi += [dj]
    mindeps[i] = depi
    #print "dep %d, %s --> %s" % (i, newdeps[i], mindeps[i])

  # 5. output the dependencies section
  sdep = []
  for i in range(n-1, 0, -1):
    if len(mindeps[i]) > 0:
      sdep += ["#ifdef %s\n" % newsrcls[i][1]]
      for j in sorted(mindeps[i], reverse = True):
         sdep += ["  #define %s\n" % newsrcls[j][1]]
      sdep += ["#endif\n", "\n"]
  return newsrcls, sdep



def mkanchors(shost, srcls):
  ''' add anchor macros to the template
  srcls is an array of (mod, MACRO) '''

  srcls, lsdeps = builddeps(srcls)

  # 1. find the location of ZCOM_PICK
  for i in range(len(shost)):
    if shost[i].startswith("#ifndef " + prefix + "PICK"):
      a0 = i + 1
      break
  else:
    print "cannot find PICK anchor"
    raise Exception

  # 2. find the location of dependence
  for i in range(a0+1, len(shost)):
    if "/* build dependencies */" in shost[i]:
      a1 = i+1
      break
  else:
    print "cannot find dependencies"
    raise Exception

  ls1 = []
  ls2 = []
  for mod, macro in srcls:
    # determine if a module is an ANSI one
    ansi = 1
    readme = os.path.join(mod, "README")
    if os.path.exists(readme):
      for s in open(readme).readlines():
        s = s.strip()
        if not s or s[0] != "#": continue
        s = s[1:].strip()
        keys = ["ADVANCED", "NONANSI", "NON-ANSI", ]
        for key in keys:
          if s.startswith(key):
            ansi = 0
            break
        if not ansi: break
    # add to ZCOM_PICK block only for an ANSI block
    if ansi:
      ls1 += ["  #ifndef %s\n" % macro, "  #define %s\n" % macro, "  #endif\n"]
    ls2 += ["#ifdef  %s\n" % macro, "#ifndef %s__\n" % macro, "#define %s__\n" %macro,
        "\n", "#endif /* %s__ */\n" % macro, "#endif /* %s */\n" % macro,
        ]  + ["\n", ] * 5

  return shost[:a0] + ls1 + shost[a0:a1] + lsdeps + shost[a1:] + ls2, srcls



def integrate(srcls, lines):
  ''' integrate a list `srcls' of source code into
      the template `lines' '''

  # 1. add anchors and dependencies
  lines, srcls = mkanchors(lines, srcls)

  modcnt = 0
  for modnm, modNM in srcls:
    if verbose > 0:
      print "add module %-8s %-13s" % (modnm, modNM)
      if verbose > 2:
        raw_input("press Enter to continue...")

    # 2. read the source code
    fnc = modnm + ".c"
    src = open(os.path.join(modnm, fnc)).readlines()
    # call rmdbg.py to remove debug/legacy code
    src = rmdbg.rmdbg(src, verbose = verbose - 1)
    # strip away the outmost #ifndef, #define, #endif triplet
    src = strip_def(src, fnc)

    # a list of files included in `src'
    for fnh in get_includes(src):
      if verbose > 1: print "%s included in %s" % (fnh, fnc)

      # 3. read the header
      hdr = open(os.path.join(modnm, fnh)).readlines()
      hdr = rmdbg.rmdbg(hdr, verbose = verbose - 1)
      hdr = strip_def(hdr, fnh)
      hdr = add_storage_class(hdr, strcls)

      # 4. insert the header to the source code
      for pivot in range(len(src)):
        if (src[pivot].startswith("#include") and
            src[pivot][8:].strip().startswith('"' + fnh + '"') ):
          break
      else:
        print "cannot find where to insert %s\n" % fnh,
        raw_input("press Enter to see the current file")
        print ''.join( src )
        raise Exception
      if verbose > 1: print "pivot is found at", pivot
      src = src[:pivot] + hdr + src[pivot + 1:]

    # 5. remove C++ (debug) commments
    src = remove_cpp_comments(src)

    # 6. insert the source code into the host
    lines = insert_module(lines, modNM, src)

    modcnt += 1

  if verbose > 0: print "%d modules, " % modcnt,

  # strip away trailing spaces
  return [ln.rstrip() + '\n' for ln in lines]



def integratef(srcls, fnhost, fnout):
  ''' integrate a list `srcls' of source code into `fnhost'
      save the result to `fnout' '''

  # 1. load the template host0
  fnr = os.path.splitext(fnhost)
  host0  = fnr[0] + ".0" + fnr[1]
  if not fnout: fnout = fnhost
  if verbose > 0:
    print "Output: %s, template: %s" % (fnout, host0)
  lines = open(host0, 'r').readlines()

  # 2. integrate
  lines = integrate(srcls, lines)

  # 3. save it back to output
  # overwrite the original only if necessary
  newsrc = ''.join(lines).strip() + "\n"
  oldsrc = ""
  if os.path.exists(fnout):
    oldsrc = open(fnout).read()
  if newsrc != oldsrc:
    print "updating", fnout
    open(fnout, 'w').write(newsrc)
  else:
    print "no need to update", fnout
  if verbose > 0:
    os.system('wc ' + fnout);



def getmodname(name):
  ''' smodule name to macro name, e.g., util --> ZCOM_UTIL '''

  shortnm = os.path.basename(name).upper()
  pos = shortnm.find('.')
  if pos >= 0: shortnm = shortnm[:pos]
  return prefix + shortnm



def getallmods(root):
  ''' get a list of modules from subdirectories of `root' '''

  # list all sub directories under defroot
  # skip modules starting with an underscore or a dot
  # also skip the test folder
  defmods = []
  for d in os.listdir(root):
    # skip regular files
    if not os.path.isdir(d): continue
    # skip folders start with special characters
    if d[0] in ('.', '_'): continue
    # skip test folder
    if d == "test": continue
    # see if the default .h and .c files exist
    file_h = os.path.join(d, d+'.h')
    if not os.path.exists(file_h):
      #print "skip %s no %s" % (d, file_h)
      continue
    file_c = os.path.join(d, d+'.c')
    if not os.path.exists(file_c):
      #print "skip %s no %s" % (d, file_c)
      continue
    defmods += [d]
  return defmods



def usage():
  """
  print usage and die
  """
  prog = sys.argv[0]
  print "%s [Options] module(s)\n" % (prog)
  print " Options:"
  print "   -a:  --all,      assemble all default modules"
  print "   -t:  --template, host file to absorb an addition"
  print "   -o:  --output,   output file"
  print "   -c:  --strcls,   storage class"
  print "   -p:  --prefix,   prefix of the host "
  print "   -v:  --verbose,  verbose, (default = 1), 0 to silence"
  print "   -h:  --help,     help\n"
  print " Example:"
  print "   %s -a" % (prog)
  print "   %s ss rng cfg -o out.h" % (prog)
  print ""
  exit(1)



def doargs():
  ''' handle arguments '''
  global fnhost, fnoutp, strcls, prefix, verbose

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hv:t:o:c:p:i:m:a",
         ["help", "prefix=", "strcls=", "template=", "output=", "src=", "all",
           "verbose=", "module="])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  defmods = getallmods(".") # util, cfg, ...

  # for each directory, guess a module name 'modNM'
  # add the pair to deflist
  deflist = []
  for fn in defmods:
    modNM = getmodname(fn) # util --> ZCOM_UTIL
    deflist += [(fn, modNM)]
  modls = []

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-a", "--all"):
      modls += deflist
    elif o in ("-t", "--template"):
      fn_host = a
    elif o in ("-o", "--output"):
      fn_outp = a
    elif o in ("-c", "--strcls"):
      strcls = a
    elif o in ("-p", "--prefix"):
      prefix = a
    elif o in ("-h", "--help"):
      usage()

  for fn in args:
    modNM = getmodname(fn)
    modls += [(fn, modNM)]

  if len(modls) == 0: usage()

  return modls



if __name__ == "__main__":
  modls = doargs() # [ (util, ZCOM_UTIL), ..., ]
  integratef(modls, fnhost, fnoutp)



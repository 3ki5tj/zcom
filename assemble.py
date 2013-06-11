#!/usr/bin/env python

'''
integrate target.c and target.h into fn_host
debug code in target.c is removed first by calling rmdbg
C++ compments are also removed
'''

import os, sys, inspect, shutil, getopt

# add the current path to the search path (unnecessary)
#curdir = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
#if curdir not in sys.path: sys.path.insert(0, curdir)

from python import rmdbg, depls


'''
module attributes, set before use
'''
fn_host     = "zcom.h"
fn_output   = None
strcls      = "STRCLS"
host_prefix = "ZCOM_"
verbose     = 1


def strip_def(src):
  '''
  strip away the
    #ifndef FILE__
    #define FILE__
    #endif
  triplet, and everything outside
  '''
  n = len(src) # number of lines
  for i in range(n):
    lin = src[i].lstrip()
    if lin.startswith("#ifndef") and (not "INLINE" in lin
        and not "RESTRICT" in lin):
      start = i
      if verbose > 1:
        print "#ifndef found in line", start
      break
  else:
    return src # no "#ifdef" found

  # look for the following "#define"
  start += 1
  if not src[start].lstrip().startswith("#define"):
    return src
  if verbose > 1:
    print "#define found in line", start

  start += 1  # beginning of the main code

  # search for the final "#endif"
  for i in range(n, start, -1):
    if src[i-1].lstrip().startswith("#endif"):
      finish = i - 1
      break
  else:
    return src

  if verbose > 1:
    print "#endif found in line", finish

  return src[start : finish]



def add_storage_class(hdr, prefix):
  '''
  add storage class, 'prefix', to global function and variables in the header
  '''
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
        "real", "float", "double", "const",
        "cfgdata_t", "logfile_t", "is_t", "hist_t",
        ) or first_word.endswith("_t"):
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



def insert_module(src, name, module):
  '''
  find the block specified by
    #ifdef  name
    #ifndef name__
    #define name__
    ...
    #endif
    #endif
  in 'src' and insert 'module' into ...
  '''
  n = len(src)
  namein = name + "__"
  plevel = 0
  bstart = -1
  bend = -1

  for i in range(n):
    if (plevel == 0 and
        bstart < 0 and
        src[i].startswith("#ifdef") and
        name in src[i] and
        src[i+1].startswith("#ifndef") and
        namein in src[i+1] and
        src[i+2].startswith("#define") and
        namein in src[i+2]):
      bstart = i+3
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

    if src[i].startswith("#if"):
      plevel += 1
    elif src[i].startswith("#endif"):
      plevel -= 1

  if bstart > 0 and bend > 0:
    return src[:bstart] + module + src[bend:]
  else:
    print "cannot find the desired module", name
    return src



def strpbrk(s, char_set):
  '''
  return a position of s where the first occurance
  of any characters in char_set is found
  '''
  try:
    pos = next(i for i,x in enumerate(s) if x in char_set)
    return pos
  except:
    return -1



def builddeps(srclist):
  ''' build dependencies '''
  n = len(srclist)
  #for i in range(n): print "%d: %s" % (i, srclist[i][0])

  # 1. build dependencies
  deps = [""]*n
  depsX = [""]*n
  for i in range(n):
    mod = srclist[i][0]
    deps[i] = []
    depsX[i] = []
    for j in range(n):
      if j == i: continue
      d = srclist[j][0]
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
      print "%d (%s) -> " % (i, srclist[i][0]),
    raise Exception

  # 3. build a new sorted src list and dependence list
  newsrclist = [0]*n
  newdeps = [0]*n
  for i in range(n):
    id = d[i]
    newsrclist[i] = srclist[id]
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
      sdep += ["#ifdef %s\n" % newsrclist[i][1]]
      for j in sorted(mindeps[i], reverse = True):
         sdep += ["  #define %s\n" % newsrclist[j][1]]
      sdep += ["#endif\n", "\n"]
  return newsrclist, sdep



def mkanchors(shost, srclist):
  ''' add anchor macros to host
  srclist is an array of (mod, MACRO) '''

  srclist, lsdeps = builddeps(srclist)

  # 1. find the location of ZCOM_PICK
  for i in range(len(shost)):
    if shost[i].startswith("#ifndef " + host_prefix + "PICK"):
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
  for mod, macro in srclist:
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
        "\n", "#endif /* %s__ */\n" % macro, "#endif /* %s */\n" % macro, "\n"]

  return shost[:a0] + ls1 + shost[a0:a1] + lsdeps + shost[a1:] + ls2, srclist



def integrate(srclist, host, fnout):
  ''' integrate source code to output according to host '''

  # 1a. load the template host0
  fnr = os.path.splitext(host)
  host0  = fnr[0] + ".0" + fnr[1]
  if not fnout: fnout = host
  if verbose > 0:
    print "Output: %s, template: %s" % (fnout, host0)
  host_src = open(host0, 'r').readlines()

  # 1b. add anchors and dependencies
  host_src, srclist = mkanchors(host_src, srclist)

  modcnt = 0

  for fn_src, mod_name in srclist:
    if not os.sep in fn_src: # make long name abc --> abc/abc
      fnlsrc = os.path.join(fn_src, fn_src)

    # derive other file names
    fnlsrc_c = fnlsrc + ".c"
    fnlsrc_h = fnlsrc + ".h"
    # short name
    fn_src     = os.path.basename(fn_src)
    fn_src_c   = fn_src + ".c"
    fn_src_h   = fn_src + ".h"
    if verbose > 1:
      print "short names are %s and %s" % (fn_src_c, fn_src_h)

    if verbose > 0:
      print ("add module %-8s %-13s to %s"
        % (fn_src, mod_name, host0))
    if verbose > 2:
      raw_input("press Enter to continue...")

    # 2. read the source code
    src = open(fnlsrc_c).readlines()
    # call rmdbg.py to remove debug/legacy code
    src = rmdbg.rmdbg(src, verbose=verbose - 1)
    # strip away the outmost #ifndef, #define, #endif triplet
    src = strip_def(src)

    # 3. read the header
    header = open(fnlsrc_h, 'r').readlines()
    header = rmdbg.rmdbg(header, verbose=verbose - 1)
    header = strip_def(header)
    header = add_storage_class(header, strcls)

    # 4. insert the header to the source code
    pivot = -1
    for i in range(len(src)):
      if (src[i].startswith("#include") and
          ('"' + fn_src_h + '"' in src[i][7:]) ):
        pivot = i
        break
    else:
      print "cannot find where to insert headers\n",
      raw_input("press Enter to see the current file")
      print ''.join( src )
      raise Exception
    if verbose > 1:
      print "pivot is found at", pivot
    src = src[:pivot] + header + src[pivot + 1:]

    # 5. remove C++ (debug) commments
    src = remove_cpp_comments(src)

    # 6. insert the source code into the host
    host_src = insert_module(host_src, mod_name, src)

    modcnt += 1
  if verbose > 0:
    print "%d modules, " % modcnt,

  for i in range(len(host_src)): # strip away trailing spaces
    host_src[i] = host_src[i].rstrip() + '\n'

  # 6. save it back to output
  # overwrite the original only if necessary
  newsrc = ''.join(host_src)
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
  if pos >= 0:
    shortnm = shortnm[:pos]
  return host_prefix + shortnm



def dirs2modules(root):
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
  global fn_host, fn_output, strcls, host_prefix, verbose

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hv:t:o:c:p:i:m:a",
         ["help", "prefix=", "strcls=", "template=", "output=", "src=", "all",
           "verbose=", "module="])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  defmods = dirs2modules(".") # util, cfg, ...

  # for each directory, guess a module name 'modnm'
  # add the pair to deflist
  deflist = []
  for fn in defmods:
    modnm = getmodname(fn) # util --> ZCOM_UTIL
    deflist += [(fn, modnm)]

  srclist = []

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-a", "--all"):
      srclist += deflist
    elif o in ("-t", "--template"):
      fn_host = a
    elif o in ("-o", "--output"):
      fn_output = a
    elif o in ("-c", "--strcls"):
      strcls = a
    elif o in ("-p", "--prefix"):
      host_prefix = a
    elif o in ("-h", "--help"):
      usage()

  for fn in args:
    modnm = getmodname(fn)
    srclist += [(fn, modnm)]

  if len(srclist) == 0:
    usage()

  return srclist

if __name__ == "__main__":
  srclist = doargs() # [ (util, ZCOM_UTIL), ..., ]
  integrate(srclist, fn_host, fn_output)



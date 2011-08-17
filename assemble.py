#!/usr/bin/env python

'''
integrate target.c and target.h into fn_host
debug code in target.c is removed first by calling rmdbg
'''

import os, shutil, getopt, sys
import rmdbg

'''
module attributes, set before use
'''
fn_host     = "zcom.h"
fn_output   = None
strcls      = "ZCSTRCLS"
host_prefix = "ZCOM_"
verbose     = 0


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
    if src[i].lstrip().startswith("#ifndef"):
      start = i
      if verbose: 
        print "#ifndef found in line", start
      break
  else:
    return src # no "#ifdef" found

  # look for the following "#define"
  start += 1
  if not src[start].lstrip().startswith("#define"):
    return src
  if verbose: 
    print "#define found in line", start

  start += 1  # beginning of the main code

  # search for the final "#endif"
  for i in range(n, start, -1):
    if src[i-1].lstrip().startswith("#endif"):
      finish = i - 1
      break
  else:
    return src
  
  if verbose: 
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
        src[i].find(name) >= 0 and 
        src[i+1].startswith("#ifndef") and
        src[i+1].find(namein) >= 0 and
        src[i+2].startswith("#define") and
        src[i+2].find(namein) >= 0):
      bstart = i+3
      if verbose:
        print "out loop starts from:", i, "line:", src[i].rstrip()

    if (plevel == 2 and
        bstart > 0 and
        src[i].startswith("#endif") and
        src[i+1].startswith("#endif") ):
      bend = i
      if verbose:
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

def checkcycl(deps):
  ''' check cyclic dependencies starting from i0 '''

  n = len(deps)
  good = [0]*n # no dependency problem
  while 1:
    # find an dep. node
    for i in range(n):
      if len(deps[i]) and not good[i]: break
    else: return 0 # all nodes are indep.

    # make a stack
    st = [-1]*(n+1)
    jt = [-1]*(n+1)
    top = 0
    st[top] = i
    jt[top] = 0 # index of child nodes of st[top]
    # tree search using a stack 
    while top >= 0:
      i = st[top]
      j = jt[top]
      dl = deps[i]
      #print "top %d: %d, deps %s, j = %d" % (top, i, dl, j); raw_input()
      if j >= len(dl): # exhausted this level
        # now anything depending on i is good
        good[i] = 1 #print "checking %d" % i
        top -= 1
        if top >= 0:
          jt[top] += 1
      elif good[ dl[j] ]:
        jt[top] += 1 # skip it
      else: # go deep into a child (dependence)
        top += 1
        x = dl[j]
        if x in st[:top]:
          print "cyclic dependence is detected from %s, stack %s" % (x, st[:top])
          return x + 1
        st[top] = x
        jt[top] = 0
        #print "adding %d on level %d" % (x, top)
  return 0

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
  i = checkcycl(depsX)
  if i != 0:
    print "Error due to cyclic dependence! at %d (%s)" % (i, srclist[i][0])
    raise Exception

  # 2. sort modules according to dependencies
  d = range(n) # ordering modules
  i = 0
  while i < n:
    #print "%d:%d(%s)  deps: %s\n%s" % (i, d[i], srclist[d[i]][0], depsX[d[i]], d)
    if len(depsX[d[i]]) == 0: 
      #print "%d:%d(%s) no more dependencies" % (i, d[i], srclist[d[i]][0]); raw_input()
      i += 1
      continue
    mi = max( map(d.index, depsX[d[i]]) )
    if mi <= i:
      i += 1
      continue
    '''
    idep = map(d.index, depsX[d[i]])
    idep.sort()
    mi = max(idep)
    if mi <= i:
      i += 1
      continue
    id = 0
    while idep[id] < i: id += 1
    mi = idep[id]
    '''
    #print "%d:%d(%s); %d:%d(%s)" % (i, d[i], srclist[d[i]][0], mi, d[mi], srclist[d[mi]][0]); raw_input()

    # swap i with mi
    d[i], d[mi] = d[mi], d[i]

  newsrclist = []
  for i in range(n):
    id = d[i]
    newsrclist += [srclist[id]]

  sdep = []
  for i in range(n-1, 0, -1):
    id = d[i]
    if len(deps[id]) > 0:
      sdep += ["#ifdef %s\n" % srclist[id][1]]
      for j in deps[id]:
         sdep += ["  #define %s\n" % srclist[j][1]]
      sdep += ["#endif\n", "\n"]
  return newsrclist, sdep

def mkanchors(shost, srclist):
  ''' add anchor macros to host '''
  
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
    if shost[i].find("build dependencies") >= 0:
      a1 = i+1
      break
  else:
    print "cannot find dependencies"
    raise Exception

  ls1 = []
  ls2 = []
  for mod, macro in srclist:
    ls1 += ["  #ifndef %s\n" % macro, "  #define %s\n" % macro, "  #endif\n"]
    ls2 += ["#ifdef  %s\n" % macro, "#ifndef %s__\n" % macro, "#define %s__\n" %macro,
        "\n", "#endif /* %s__ */\n" % macro, "#endif /* %s */\n" % macro, "\n"]

  return shost[:a0] + ls1 + shost[a0:a1] + lsdeps + shost[a1:] + ls2, srclist


def integrate(srclist, host, output):
  ''' integrate source code to output according to host '''

  # 1. load the template host0
  fnr = os.path.splitext(host)
  host0  = fnr[0] + ".0" + fnr[1]
  if not output: output = host
  print "Output: %s, template: %s" % (output, host0)
  host_src = open(host0, 'r').readlines()
  host_src, srclist = mkanchors(host_src, srclist)

  modcnt = 0

  for fn_src, mod_name in srclist:
    if fn_src.find(os.sep) < 0: # make long name abc --> abc/abc
      fnlsrc = os.path.join(fn_src, fn_src)
  
    # derive other file names
    fnlsrc_c = fnlsrc + ".c"
    fnlsrc_h = fnlsrc + ".h"
    # short name
    fn_src     = os.path.basename(fn_src)
    fn_src_c   = fn_src + ".c"
    fn_src_h   = fn_src + ".h"
    if verbose:
      print "short names are %s and %s" % (fn_src_c, fn_src_h)
  
    print ("add module %-8s %-13s to %s" 
        % (fn_src, mod_name, host0))
    if verbose:
      raw_input("press Enter to continue...")
 
    # 2. read the source code
    src = open(fnlsrc_c).readlines()
    # call rmdbg.py to remove debug/legacy code
    src = rmdbg.rmdbg(src, verbose=verbose)
    # strip away the outmost #ifndef, #define, #endif triplet
    src = strip_def(src)
  
    # 3. read the header
    header = open(fnlsrc_h, 'r').readlines()
    header = rmdbg.rmdbg(header, verbose=verbose)
    header = strip_def(header)
    header = add_storage_class(header, strcls)
  
    # 4. insert the header to the source code
    pivot = -1
    for i in range(len(src)):
      if (src[i].startswith("#include") and 
          src[i].find('"' + fn_src_h + '"', 7) >= 0):
        pivot = i
        break
    else:
      print "cannot find where to insert headers\n", 
      raw_input("press Enter to see the current file")
      print ''.join( src )
      return 1
    if verbose: 
      print "pivot is found at", pivot
    src = src[:pivot] + header + src[pivot + 1:]
    
    # 5. insert the source code into the host
    host_src = insert_module(host_src, mod_name, src)

    modcnt += 1
  print "%d modules, " % modcnt,

  for i in range(len(host_src)): # strip away trailing spaces
    host_src[i] = host_src[i].rstrip() + '\n'
  
  # 6. save it back to output
  # save first to a temporary file,
  # overwrite the original if necessary
  fn_outtmp = output + ".tmp"
  open(fn_outtmp, 'w').write(''.join(host_src))
  if os.system('cmp '+fn_outtmp + ' '+ output):
    shutil.copy(fn_outtmp, output)
  else:
    print "no need to update", output
  os.remove(fn_outtmp)

  os.system('wc ' + output);

def getmodname(name):
  ''' smodule name to macro name, e.g., util --> ZCOM_UTIL '''
  shortnm = os.path.basename(name).upper()
  pos = shortnm.find('.')
  if pos >= 0:
    shortnm = shortnm[:pos]
  return host_prefix + shortnm

def dirs2mods(root):
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
    # see if a default C module exists
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
  print "   -a:  --all,      all default modules"
  print "   -t:  --template, host file to absorb an addition"
  print "   -o:  --output,   output file"
  print "   -c:  --strcls,   storage class"
  print "   -p:  --prefix,   prefix of the host "
  print "   -v:  --verbose,  verbose"
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
 
  defmods = dirs2mods(".") # util, cfg, ...
  
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
  print "modules:", srclist
  integrate(srclist, fn_host, fn_output)



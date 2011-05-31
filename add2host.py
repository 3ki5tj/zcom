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
        "real", "float", "double", "cfgdata_t", "logfile_t", "is_t"):
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

def get_mod_name(name):
  '''
  set the module name inside host
  '''
  global mod_name

  if name.startswith(host_prefix):
    modname = name
  else: 
    # make a guess
    shortnm = os.path.basename(name).upper()
    pos = shortnm.find('.')
    if pos >= 0:
      shortnm = shortnm[:pos]
    modname = host_prefix + shortnm
  return modname


def integrate(srclist):
  '''
  integrate fn_source to fn_host (output)
  according to a template fn_host_t (with .0 extension)
  '''

  # 1. load the template fn_host_t 
  fnr = os.path.splitext(fn_host)
  fn_host_t  = fnr[0] + ".0" + fnr[1]
  print "Host: %s, template: %s" % (fn_host, fn_host_t)
  host_src = open(fn_host_t, 'r').readlines()

  for fn_source, mod_name in srclist:
    if fn_source.find(os.sep) < 0:
      # make abc --> abc/abc
      fn_source = os.path.join(fn_source, fn_source)
  
    # derive other file names
    fn_source_c = fn_source + ".c"
    fn_source_h = fn_source + ".h"
    fn_host_bak = fn_host + ".bak"
    # short name
    fn_src     = os.path.basename(fn_source)
    fn_src_c   = fn_src + ".c"
    fn_src_h   = fn_src + ".h"
    if verbose:
      print "short names are %s and %s" % (fn_src_c, fn_src_h)
  
    print ("integrating module %-12s (%-18s, %-18s) to %s" 
        % (mod_name, fn_source_c, fn_source_h, fn_host_t))
    if verbose:
      raw_input("press Enter to continue...")
  
  
    # 2. read the source code
    src = open(fn_source_c, 'r').readlines()
    # call rmdbg.py to remove debug/legacy code
    src = rmdbg.rmdbg(src, verbose=verbose)
    # strip away the outmost #ifndef, #define, #endif triplet
    src = strip_def(src)
  
    # 3. read the header
    header = open(fn_source_h, 'r').readlines()
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

  # 6. save it back to fn_host
  # save first to a temporary file,
  # overwrite the original if necessary
  fn_host_tmp = fn_host + ".tmp"
  open(fn_host_tmp, 'w').write(''.join(host_src))
  if os.system('cmp '+fn_host_tmp + ' '+ fn_host):
    shutil.copy(fn_host_tmp, fn_host)
  else:
    print "no need to update", fn_host
  os.remove(fn_host_tmp)

  os.system('wc ' + fn_host);

def get_defmodules(root):
  '''
  obtain a list of C modules under `root'
  '''
  # list all sub directories under defroot
  # skip modules starting with an underscore or a dot
  # also skip the test folder
  defmodules = []
  for d in os.listdir(root):
    # skip regular files
    if not os.path.isdir(d): continue
    # skip folders start with special characters
    if d[0:1] in ('.', '_'): continue
    # skip test folder
    if d == "test": continue
    # see if a default C module exists
    file_h = os.path.join(d, d+'.h')
    if not os.path.exists(file_h): 
      print "skip", d, "no", file_h
      continue
    file_c = os.path.join(d, d+'.c')
    if not os.path.exists(file_c): 
      print "skip", d, "no", file_c
      continue
    defmodules += [d]
  return defmodules


def usage():
  """
  print usage and die
  """
  print sys.argv[0], "[Options] module(s)"
  print " Options:"
  print " -a:  --all,     all default modules"
  print " -o:  --host,    host file to absorb an addition"
  print " -c:  --strcls,  storage class"
  print " -p:  --prefix,  prefix of the host "
  print " -v:  --verbose, verbose\n"
  print " Example:"
  print sys.argv[0], "ss rng cfg"
  exit(1)

def handle_params():
  '''
  Handle common parameters from command line options
  results saved to module attributes
  '''
  global fn_host, fn_source, strcls, host_prefix, verbose, mod_name

  try:
    opts, args = getopt.getopt(sys.argv[1:], "h:v:o:c:p:i:m:a", 
         ["help", "prefix=", "strcls=", "host=", "src=", "all", 
           "verbose=", "module="])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
  
  defmodules = get_defmodules(".")
  
  # for each directory, guess a module name 'modnm'
  # add the pair to deflist
  deflist = []
  for fn in defmodules: 
    modnm = get_mod_name(fn)
    deflist += [(fn, modnm)]
  
  srclist = []

  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-a", "--all"):
      srclist += deflist
    elif o in ("-h", "--host", "-o"):
      fn_host = a
    elif o in ("-c", "--strcls"):
      strcls = a
    elif o in ("-p", "--prefix"):
      host_prefix = a
    elif o in ("-m", "--module"):
      mod_name = a
    elif o in ("-h", "--help"):
      usage()

  for fn in args:
    modnm = get_mod_name(fn)
    srclist += [(fn, modnm)]

  if len(srclist) == 0:
    usage()

  return srclist

def main():
  '''
  driver
  '''
  srclist = handle_params()
  print "modules:", srclist
  integrate(srclist)

if __name__ == "__main__":
  main()



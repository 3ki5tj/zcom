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
verbose     = 1

fn_source   = ""    # e.g., "ss/ss"
mod_name    = ""    # name in host

def strip_def(src, verbose = 1):
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
    arr = hdr[i].split()
    if len(arr) > 0:
      first_word = arr[0]
    else:
      continue
    
    #  NOTE: we cannot handle more complex cases
    if first_word in ("void", "char", "int", "unsigned", "long", "float", "double"):
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

def set_mod_name(name):
  '''
  set the module name inside host
  '''
  global mod_name

  if name.startswith(host_prefix):
    mod_name = name
  else: 
    # make a guess
    shortnm = os.path.basename(name).upper()
    pos = shortnm.find('.')
    if pos >= 0:
      shortnm = shortnm[:pos]
    mod_name = host_prefix + shortnm
  return mod_name


def integrate():
  '''
  integrate fn_source to fn_host
  '''
  
  # make a local copy
  fn_src_l = fn_source

  if fn_src_l.find(os.sep) < 0:
    # make abc --> abc/abc
    fn_src_l = os.path.join(fn_src_l, fn_src_l)

  # derive other file names
  fn_source_c = fn_src_l + ".c"
  fn_source_h = fn_src_l + ".h"
  fn_host_bak = fn_host + ".bak"
  # short name
  fn_src     = os.path.basename(fn_src_l)
  fn_src_c   = fn_src + ".c"
  fn_src_h   = fn_src + ".h"
  if verbose:
    print "short names are %s and %s" % (fn_src_c, fn_src_h)

  print ("integrating module %s (%s, %s) to %s" 
      % (mod_name, fn_source_c, fn_source_h, fn_host))
  if verbose:
    raw_input("press Enter to continue...")


  # 1. read the source code
  src = open(fn_source_c, 'r').readlines()
  # call rmdbg.py to remove debug information
  src = rmdbg.rmdbg(src)
  # strip away the outmost #ifndef, #define, #endif triplet
  src = strip_def(src)

  # 2. read the header
  header = open(fn_source_h, 'r').readlines()
  header = strip_def(header)
  header = add_storage_class(header, strcls)

  # 3. insert the header to the source code
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

  # 4. load fn_host and insert the source into it 
  host_src = open(fn_host, 'r').readlines()
  # insertion
  host_src = insert_module(host_src, mod_name, src)
  
  # 5. save it back to fn_host
  shutil.copy2(fn_host, fn_host_bak) # make a backup first
  open(fn_host, 'w').write(''.join(host_src))

def usage():
  """
  print usage and die
  """
  print sys.argv[0], "[Options]"
  print "Options:"
  print " -s:  --src,     source code to be integrated"
  print " -o:  --host,    host file to absorb an addition"
  print " -c:  --strcls,  storage class"
  print " -p:  --prefix,  prefix of the host "
  print " -v:  --verbose, verbose"
  exit(1)

def handle_params():
  '''
  Handle common parameters from command line options
  results saved to module attributes
  '''
  global fn_host, fn_source, strcls, host_prefix, verbose, mod_name

  try:
    opts, args = getopt.getopt(sys.argv[1:], "h:v:s:o:c:p:i:m:", 
         ["help", "prefix=", "strcls=", "host=", "src=", 
           "verbose=", "module="])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
  
  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = int(a)
    elif o in ("-s", "--src", "-i"):
      fn_source = a
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
 
  if fn_source == "":
    usage()

  if mod_name == "":
    set_mod_name(fn_source)

def main():
  '''
  driver
  '''
  handle_params()
  integrate()

if __name__ == "__main__":
  main()


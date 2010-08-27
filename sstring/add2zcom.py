#!/usr/bin/env python

'''
integrate ss.c and ss.h into zcom.c
currently assume the debug information has already been stripped away
'''

import os, shutil;

def strip_def(src, verbose=1):
  '''
  strip away the #ifndef, #define, #endif triplet, and things outside them
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


def insert_module(src, name, module, verbose=1):
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
    print "cannot find the desired module"
    return src

def main(verbose=1):
  fn_ss_c = "ss.c"
  fn_zcom_c = "../zcom.c"
  strcls = "ZCSTRCLS"

  # derive other file names
  fn_ss_h = os.path.splitext(fn_ss_c)[0] + ".h"
  fn_zcom_bak = fn_zcom_c + ".bak"

  # call rmdbg.py to remove debug information
  # ...

  # 1. read the source code
  f=open(fn_ss_c, 'r')
  src=f.readlines()
  f.close()
  # strip away the #ifndef, #define, #endif triplet
  src = strip_def(src)
  # print ''.join(src)

  # 2. read the header
  f=open(fn_ss_h, 'r')
  header = f.readlines()
  f.close()
  header = add_storage_class(strip_def(header), strcls)

  # 3. insert the header to the source code
  pivot = -1
  for i in range(len(src)):
    if (src[i].startswith("#include") and 
        src[i].find('"ss.h"', 7) >= 0):
      pivot = i
      break
  if verbose:
    if pivot < 0:
      print "cannot find where to insert headers"
    else:
      print "pivot is found at", pivot
  src = src[:pivot] + header + src[pivot+1:]
  #print ''.join(src)

  # 4. load zcom.c and insert the source into it 
  f=open(fn_zcom_c, 'r')
  zcom = f.readlines()
  f.close()
  # insertion
  zcom = insert_module(zcom, "ZCOM_SS", src)
  #print ''.join(zcom)
  
  # 5. save it back to zcom.c
  # make a backup
  shutil.copy2(fn_zcom_c, fn_zcom_bak)
  f = open(fn_zcom_c, 'w')
  f.write(''.join(zcom))
  f.close()

if __name__ == "__main__":
    main()


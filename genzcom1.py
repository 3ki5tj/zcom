#!/usr/bin/env python

# try to create a limited version of zcom1.h from zcom.h

# read the file 
f=open('zcom.h', 'r')
in_lines=f.readlines()
f.close()

verbose=0

# briefly search all keys in the toc section
# also exclude all good keys
def find_keys(lines, excls):
  start_toc=-1
  depth=0
  list=[]
  lnum=1
  for line in lines:
    lin1=line.strip()
    if lin1.startswith("#if"):
      depth+=1

    if lin1.startswith("#ifndef"):
      if depth==1 and (start_toc<0) and lin1.find("ZCOM_PICK")>=0:
        start_toc=depth
      elif start_toc>=0 and depth==2:
        newkey=lin1.split()[1]
        if not (newkey in excls):
          list += [ newkey ]

    if lin1.startswith("#endif"):
      if start_toc>=0 and start_toc==depth:
        start_toc=-1
        return list
      depth-=1
    
    lnum+=1

  return list

good_keys=['ZCOM_RNG','ZCOM_CFG','ZCOM_TRACE','ZCOM_DIST','ZCOM_ENDIAN','ZCOM_ERROR']
bad_keys=find_keys(in_lines, good_keys)
print "good:", good_keys
print "bad :", bad_keys


out_lines=[]
plevel=0
linenum=1
lnumout=1

ignore_at=-1
toc_at=-1
inside_readme=0

last_nonblank=0

def has_key(line, keys):
  for key in keys:
    if line.find(key) >= 0:
      if verbose >= 2:
        print "IGNOR +",plevel,"linenum=",linenum, ":", line.rstrip(), "key:", key
      return 1
  return 0

for line in in_lines:
  lin1=line.lstrip()
  
  # readme is the first block comment in the file
  if inside_readme==0:
    if lin1.startswith("/*"):
      inside_readme=1
  elif inside_readme>0:
    # the splitter for additional informations
    if inside_readme==1 and lin1.startswith("-------"):
      inside_readme=2
    elif lin1.find("*/")>=0:
      inside_readme=-1

  # increase the level
  if lin1.startswith("#if"):
    plevel+=1
    if verbose>=3:
      print "LEVEL +",plevel,"linenum=",linenum, ":", lin1

    # the "table of contents" led by ZCOM_PICK
    if toc_at<0:
      if lin1.startswith("#ifndef ZCOM_PICK"):
        toc_at=plevel
    elif ignore_at<0 and lin1.startswith("#ifndef"):
      if has_key(lin1, bad_keys): ignore_at=plevel

    # ignore
    if ignore_at<0:
      if plevel==1 and lin1.startswith("#ifdef"):
        if has_key(lin1, bad_keys): ignore_at=plevel
      elif plevel==1 and lin1.find("ZCOM_RVDIM_BINDING")>=0:
        ignore_at=plevel
      elif lin1.find("_LEGACY")>=0:
        ignore_at=plevel
        
  if verbose>=3:
    print "#", linenum, "plevel:",plevel,"ignore_at:",ignore_at,"R:",inside_readme,"line:", lin1.rstrip()
  

  if ignore_at<0 and inside_readme<2:
    if lin1!="": last_nonblank=lnumout

    # we allow two successive blank lines, but no more
    if lnumout-last_nonblank<=2:
      out_lines += [line]
      lnumout+=1
  

  # we update plevel change due #endif here
  if lin1.startswith("#endif"):
    if toc_at >= 0 and toc_at == plevel:
      toc_at=-1

    if ignore_at >= 0 and ignore_at == plevel:
      if verbose >= 2:
        print "IGNOR -",plevel,"linenum=",linenum, ":", lin1
   
      ignore_at=-1

    plevel-=1
    if verbose>=3:
      print "LEVEL -",plevel,"linenum=",linenum, ":", lin1  
    if(plevel<0):
      print "plevel is negative at line:", linenum
      exit(1)
  
  # increase the line number
  linenum+=1
  

# output file
f=open('zcom1.h', 'w')
f.write(''.join(out_lines))
f.close()




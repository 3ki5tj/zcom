#!/usr/bin/env python

# try to create a limited version of zcom1.h from zcom.h

# read the file 
f=open('ss.c', 'r')
in_lines=f.readlines()
f.close()

verbose=3

bad_keys=['SSDBG_']
print "bad :", bad_keys

out_lines=[]
plevel=0
linenum=1
lnumout=1
ignore_at=-1
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
  
  # increase the level
  if lin1.startswith("#if"):
    plevel+=1
    if verbose>=3:
      print "LEVEL +",plevel,"linenum=",linenum, ":", lin1

    # ignore
    if ignore_at<0:
      if lin1.startswith("#ifdef"):
        if has_key(lin1, bad_keys): ignore_at=plevel
        
  if verbose>=3:
    print "#", linenum, "plevel:",plevel,"ignore_at:",ignore_at,lin1.rstrip()
  

  if ignore_at<0:
    if lin1!="": last_nonblank=lnumout

    # we allow two successive blank lines, but no more
    if lnumout-last_nonblank<=2:
      out_lines += [line]
      lnumout+=1
  

  # we update plevel change due #endif here
  if lin1.startswith("#endif"):
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
f=open('ss.clean.c', 'w')
f.write(''.join(out_lines))
f.close()




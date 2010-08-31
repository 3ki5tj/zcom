#!/usr/bin/env python

"""
remove debug code from ss.c
"""

def rmdbg(in_lines, badkeys=['SSDBG_', 'CFGDBG_', '_LEGACY'], verbose=2):
  """
  read lines from instr, remove badkeys, 
  return 
  """
  print "The debug code signatures are", badkeys
  
  out_lines=[]
  plevel=0
  linenum=1
  lnumout=1
  ignore_at=-1
  last_nonblank=0

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
          # try to scan over badkeys
          for key in badkeys:
            if lin1.find(key) >= 0:
              ignore_at = plevel
              if verbose >= 2:
                print "IGNOR +",plevel,"linenum=",linenum, ":", line.rstrip(), "key:", key
              break
          else:
            pass # do nothing

    if verbose>=3:
      print "#", linenum, "plevel:", plevel, "ignore_at:",ignore_at,lin1.rstrip()
  
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
  
  return out_lines

def main(input, output):
  """
  read file 'input' 
  remove debug parts, as signified by badkeys,
  save result to 'output'
  also remove three or more successive blank lines
  """

  # read input, call rmdbg, and write to output
  open(output, 'w').write(''.join( rmdbg(open(input, 'r').readlines()) ))

if __name__ == "__main__":
  main('ss.c', 'ss.clean.c')


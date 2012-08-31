#!/usr/bin/env/python
import os, sys, runcom

def sizedep_mcl(fnrep):
  ''' Monte Carlo, local for the size dependence '''
  runcom.todatadir()
  fnlog = "sizedep_mcl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd0 = "..\ljmc -w -1 100000" # common parameters
  lines = ["# perturb-size      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  for sz in (0.001, 0.002, 0.005, 0.01, 0.02):
    cmd = cmd0
    # set size and histogram range
    cmd += " -M %s" % sz # set size
    cmd += " -u %s -U %s " % (sz*.5, sz*5000) # set histogram
    cmd += " -o ehmclM%s.dat" % sz # output
    cmd += " -O ehmcldM%s.dat" % sz # output
    d, lns = runcom.getoutp(cmd, fnlog)
    s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  0" % (sz,
        d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
        d["bph0"], d["bph1"], d["bpd0"], d["bpd1"], d["bpdi"],
        d["bc"], d["bc0"], d["bcr"])
    lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")

sizedep_mcl("sizemcl.txt")


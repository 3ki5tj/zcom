#!/usr/bin/env python
import os, sys, runcom

nsteps = 100000

def sizedep_mcl(fnrep):
  ''' Monte Carlo, local for the size dependence '''
  runcom.todatadir()
  fnlog = "sizedep_mcl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd0 =  os.path.join("..", "ljmc") + " -w -1 %s" % nsteps
  lines = ["# perturb-size      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  for sz in (0.005, 0.01, 0.02, 0.05, 0.1):
    cmd = cmd0
    # set size and histogram range
    cmd += " -M %s" % sz # set size
    cmd += " -u %s -U %s " % (sz*.1, sz*1000) # set histogram
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

if len(sys.argv) > 1:
  nsteps = int(sys.argv[1])
print "simulation length: %s" % nsteps

sizedep_mcl("sizemcl.txt")


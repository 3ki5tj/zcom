#!/usr/bin/env python
import os, sys, runcom

nsteps = 100000

def rcdep_mcl(fnrep):
  ''' Monte Carlo, local for the cutoff dependence '''
  runcom.todatadir()
  fnlog = "rcdep_mcl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd0 = os.path.join("..", "ljmc") + " -1 %s" % nsteps
  sz = 0.05
  cmd0 += " -M %s" % sz # set size
  cmd0 += " -u %s -U %s " % (sz*.1, sz*1000) # set histogram  
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  for rc in (1.2, 1.4, 1.6, 2.0, 2.5):
    cmd = cmd0
    # set cutoff
    cmd += " -c %s -s %s" % (rc, rc) # set cutoff
    cmd += " -o ehmclc%s.dat" % rc # output
    d, lns = runcom.getoutp(cmd, fnlog)
    s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  0" % (rc,
        d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
        d["bph0"], d["bph1"], d["bpd0"], d["bpd1"], d["bpdi"],
        d["bc"], d["bc0"], d["bcr"])
    lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")

if len(sys.argv) > 1:
  nsteps = int(sys.argv[1])
print "simulation length: %s" % nsteps

rcdep_mcl("rcmcl.txt")


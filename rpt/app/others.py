#!/usr/bin/env python
import os, sys, runcom

nsteps = 100000

def mcg(fnrep):
  ''' Monte Carlo, global move '''
  runcom.todatadir()
  fnlog = "mcg.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd = os.path.join("..", "ljmc") + " -w -g -1 %s" % nsteps
  sz = 0.01
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  cmd += " -M %s" % sz # set size
  cmd += " -u %s -U %s " % (sz*1, sz*10000) # set histogram  
  cmd += " -o ehmcgM%s.dat" % sz # output
  cmd += " -O ehmcgdM%s.dat" % sz # output
  d, lns = runcom.getoutp(cmd, fnlog)
  s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  1" % (sz,
      d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
      d["bpd0"], d["bpd1"], d["bpdi"],
      d["bc"], d["bc0"], d["bcr"])
  lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")

def mdl(fnrep):
  ''' Molecular dynamics, local move '''
  runcom.todatadir()
  fnlog = "mdl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd = os.path.join("..", "ljmd") + " -w -1 %s" % nsteps
  sz = 0.05
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  cmd += " -M %s" % sz # set size
  cmd += " -u %s -U %s " % (sz*.1, sz*1000) # set histogram  
  cmd += " -o ehmdlM%s.dat" % sz # output
  cmd += " -O ehmdldM%s.dat" % sz # output
  d, lns = runcom.getoutp(cmd, fnlog)
  s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  0" % (sz,
      d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
      d["bph0"], d["bph1"], d["bpd0"], d["bpd1"], d["bpdi"],
      d["bc"], d["bc0"], d["bcr"])
  lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")


def sqmcl(fnrep):
  ''' square-well, Monte Carlo, local move '''
  runcom.todatadir()
  fnlog = "sqmcl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd0 = os.path.join("..", "sqmc") + " -1 %s" % nsteps
  for sz in (0.05, 0.1, 0.2):
    lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
    cmd = cmd0
    cmd += " -M %s" % sz # set size
    cmd += " -o ehsqmclM%s.dat" % sz # output
    cmd += " -O ehsqmcldM%s.dat" % sz # output
    d, lns = runcom.getoutp(cmd, fnlog)
    s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  0" % (sz,
        d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
        d["bph0"], d["bph1"], d["bpd0"], d["bpd1"], d["bpdi"])
    lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")


def ismcl(fnrep):
  ''' Ising, Monte Carlo, local move '''
  runcom.todatadir()
  fnlog = "ismcl.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd = os.path.join("..", "is") + " -1 %s" % nsteps
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  cmd += " -o ehismcl.dat" # output
  cmd += " -O ehismcld.dat" # output
  d, lns = runcom.getoutp(cmd, fnlog)
  s = "%-18s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  0" % (4,
      d["bp0"], d["bp1"], d["bpi"], d["bps0"], d["bps1"],
      d["bph0"], d["bph1"], d["bpd0"], d["bpd1"], d["bpdi"])
  lines += [ s ]
  open(fnrep, "w").write('\n'.join(lines) + '\n')
  os.chdir("..")

def isent():
  fnlog = "data/isent.log"
  open(fnlog, "w").write("") # clear the log file
  size = 32
  fnprof = "data/profis%s.dat" % size
  fnflow = "data/flowis%s.dat" % size
  cmd = os.path.join(".", "isent") + " -o %s -f %s -1 %s" % (fnprof, fnflow, nsteps)
  d, lns = runcom.getoutp(cmd, fnlog)

if len(sys.argv) > 1:
  nsteps = int(sys.argv[1])
print "simulation length: %s" % nsteps


mcg("mcg.txt")

mdl("mdl.txt")

sqmcl("sqmcl.txt")

#ismcl("ismcl.txt")


#isent()


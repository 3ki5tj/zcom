#!/usr/bin/env python
import os, sys, runcom

nsteps = 1000000

def mcg(fnrep):
  ''' Monte Carlo, global move '''
  runcom.todatadir()
  fnlog = "mcg.log"
  open(fnlog, "w").write("") # clear the log file 
  cmd = "..\ljmc -w -g -1 %s" % nsteps
  sz = 0.002
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  cmd += " -M %s" % sz # set size
  cmd += " -u %s -U %s " % (sz*5, sz*50000) # set histogram  
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
  cmd = "..\ljmd -w -1 %s" % nsteps
  sz = 0.01
  lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
  cmd += " -M %s" % sz # set size
  cmd += " -u %s -U %s " % (sz*.5, sz*5000) # set histogram  
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
  cmd = "..\sqmc -1 %s" % nsteps
  for sz in (0.01, 0.02, 0.04):
    lines = ["# cutoff-dist.      bp0       bp1       bpi       bps0      bps1      "
      "bph0      bph1      bpd0      bpd1      bpdi      bc        bc0       bcr       g",]
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
  cmd = "..\is -1 %s" % nsteps
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
  fnprof = "profis%s.dat" % size
  fnflow = "flowis%s.dat" % size
  cmd = "..\isent -o %s -f %s" % (fnprof, fnflow)
  d, lns = runcom.getoutp(cmd, fnlog)
  shutil.move(fnprof, os.path.join("data", fnprof))
  shutil.move(fnflow, os.path.join("data", fnflow))

mcg("mcg.txt")

mdl("mdl.txt")

sqmcl("sqmcl.txt")

ismcl("ismcl.txt")

#isent()


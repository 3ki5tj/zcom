#!/usr/bin/env python
''' collects data for the Lennard Jones system '''
import os, sys

def getoutp(cmd, fnlog = None):
  #tmpfn = "TMPOUTPT.LOG"
  tmpfn = os.tempnam()
  print "CMD:", cmd
  os.system("%s > %s" % (cmd, tmpfn))
  if not os.path.exists(tmpfn):
    print "no output for %s" % cmd
    return None, None
  
  lines = open(tmpfn).read().strip().split("\n")
  lines = [cmd] + lines

  # parse the last line to a dictionary, which looks like
  # a 1, bc 3.2, efg 4.5
  data = lines[-1].strip().split(",")
  d = {}
  for dat in data:
    try:
      kd = dat.strip().split()
      d[kd[0]] = kd[1]
    except Exception:
      print "error", dat, "|", kd

  if fnlog:
    open(fnlog, "a").write('\n'.join(lines) + "\n\n")
  os.remove(tmpfn)
  return d, lines


def todatadir(dir = "data"):
  ''' change the current directory to `data' '''
  if not os.path.exists(dir):
    os.mkdir(dir)
  os.chdir(dir)

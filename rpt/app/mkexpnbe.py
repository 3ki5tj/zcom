#!/usr/bin/env python
import os, sys, glob
from math import *

def getexp(his, bet):
  cnt = 0
  smxp = 0
  for k in his:
    h = his[k]
    cnt += h
    xp = exp(-bet*k)
    smxp += xp * h
  return smxp / cnt

def mkexpnbe(fn):
  ''' beta vs. exp(-beta * e) '''
  lines = open(fn).read().strip().split('\n')
  if lines[0].startswith("#"):
    info = lines[0].split()
    lines = lines[1:]
  his = {} # using dict to construct a histogram
  for l in lines:
    a = l.split()
    #print a
    if a[0].startswith("inf"): continue # infinity
    e = float(a[0])
    h = float(a[1])
    if h <= 0.: continue
    his[e] = h

  m = 100
  dbet = 1.0/m
  s = ""
  for ib in range(-m, 2*m+1): 
    bet = ib*dbet
    s += "%.6f %.6f\n" % (bet, getexp(his, bet))

  fnout = os.path.splitext(fn)[0] + ".eb"
  open(fnout, "w").write(s)
  print "saved output to", fnout



if len(sys.argv) <= 1:
  fnls = [ "eh*.dat" ]
else:
  fnls = sys.argv[1:]

ls = []
for ar in fnls:
  ls += glob.glob(ar)

for fn in ls:
  mkexpnbe(fn)



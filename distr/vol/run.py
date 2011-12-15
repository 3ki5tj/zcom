#!/usr/bin/env python

''' run the program eh with different conditions, output data '''

import os, sys, re, shutil

active = 0 # run md to collect data
datadir = "data"
fngp = "vol.gp"
fncfg = "vol.cfg"

def foo(runmd = 0):
  # write configure file
  cfgstr = open(os.path.join("..", fncfg)).readlines()[2:]
  if runmd:
    cfgstr = ["initload = 0\n", "dosimul = 1\n"] + cfgstr
  else:
    cfgstr = ["initload = 1\n", "dosimul = 0\n"] + cfgstr
  open(fncfg, "w").writelines(cfgstr)

  # run vol
  os.system("../vol")

''' Main program starts here '''
# build the c program
os.system("make")

# change data dir
os.chdir(datadir)

foo(runmd = active)

# call Gnuplot to do the ploting
try:
  import Gnuplot
  g = Gnuplot.Gnuplot(persist = 1)
  g.load(fngp)
  g.reset()
except: # no gnuplot module
  print "cannot find gnuplot module for python, skip plotting"


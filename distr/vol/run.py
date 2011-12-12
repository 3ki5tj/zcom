#!/usr/bin/env python

''' run the program eh with different conditions, output data '''

import os, sys, re, shutil

active = 1
datadir = "data"
fngp = "vol.gp"

def run(run = 0):
  pass

''' Main program starts here '''
# build the c program
os.system("make")

# change data dir
os.chdir(datadir)

run(run = active)

# call Gnuplot to do the ploting
try:
  import Gnuplot
  g = Gnuplot.Gnuplot(persist = 1)
  g.load(fngp)
  g.reset()
except: # no gnuplot module
  print "cannot find gnuplot module for python, skip plotting"


#!/usr/bin/env python

''' run the program eh with different conditions, output data '''

import os, sys, re, shutil

active = 0 # run md to collect data
datadir = "data"
fngp = "rdf.gp"
fngpref = "rdfref.gp"
fncfg = "rdf.cfg"

def foo(runmd = 0, tp = 0.85, aj = 0):
  print "tp %g aj %d" % (tp, aj)

  # write configure file
  cfgstr = open(os.path.join("..", fncfg)).readlines()[2:]
  '''
  # use the following code to change other settings
  for k in range(len(cfgstr)):
    if cfgstr[k].startswith("nsteps"):
      cfgstr[k] = "nsteps = 1000000\n"
    elif cfgstr[k].startswith("nstdb"):
      cfgstr[k] = "nstdb = 100\n"
  '''

  if runmd:
    cfgstr = ["initload = 0\n", "dosimul = 1\n"] + cfgstr
  else:
    shutil.copy2("ds%s.dat" % tp, "ds.dat")
    shutil.copy2("dsb%s.dat" % tp, "dsb.dat")
    cfgstr = ["initload = 1\n", "dosimul = 0\n"] + cfgstr

  for i in range(len(cfgstr)):
    if cfgstr[i].startswith("iitype"):
      if aj:
        cfgstr[i] = "iitype = 9\n"
      else:
        cfgstr[i] = "iitype = 1\n"
    elif cfgstr[i].startswith("tp"):
      cfgstr[i] = "tp = %s\n" % tp
    
  open(fncfg, "w").writelines(cfgstr)

  # run rdf
  os.system("../rdf")

  # file name tag
  tag = ""
  if aj: tag = "aj"
  tag += "%s" % tp
  
  # copy
  shutil.copy2("ds.dat", "ds%s.dat" % tag)
  shutil.copy2("dsb.dat", "dsb%s.dat" % tag)
  shutil.copy2("rdf.cfg", "rdf%s.cfg" % tag)

''' Main program starts here '''
# build the c program
os.system("make")

# change data dir
if not os.path.isdir(datadir):
  os.mkdir(datadir)
os.chdir(datadir)
os.system("ln -s ../*.gp .")

foo(runmd = active, tp = 0.85)
foo(runmd = 0, tp = 0.85, aj = 1)
foo(runmd = active, tp = 0.4)
foo(runmd = 0, tp = 0.4, aj = 1)

# call Gnuplot to do the ploting
if fngp:
  try:
    import Gnuplot
    g = Gnuplot.Gnuplot(persist = 1)
    g.load(fngp)
    g.reset()
    g.load(fngpref)
  except: # no gnuplot module
    print "cannot find gnuplot module for python, skip plotting"


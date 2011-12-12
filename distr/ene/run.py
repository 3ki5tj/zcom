#!/usr/bin/env python

''' run the program eh with different conditions, output data '''

import os, sys, re, shutil

active = 1
datadir = "data"
fngp = "ene.gp"

# template for eh_cfg
src_cfg = '''initload = %s
dosimul = %s
nequil = 20000
nstadj = 2000
nsteps = 1000000
nstdb = 100 # nst of saving to the smaller database dsb
tstat = %s
n = 256
rho = 0.8
rc = 3.0
rs = 2.0
tp = 1.0
emin = -1400.0
emax = -1220.0
edel = 0.1
iitype = %-4s     # 0: AJ identity; 1: fractional identity
halfwin = %-4s    # 0: single bin; > 0: fixed win; -1: variable window
mfhalfwin = %-4s  # 0: signle bin; > 0: ii; < 0: plain average
'''

def run(canon = 1, run = 0, iitype = 1, halfwin = 70, mfhalfwin = 0):
  # write a config file
  if run:
    initload = 0
    dosimul = 1
  else:
    initload = 1
    dosimul = 0
  strcfg = src_cfg % (initload, dosimul, canon, iitype, halfwin, mfhalfwin)

  # make file names
  if canon:
    fnds = "ds_c"
    fndsb = "dsb_c"
    fncfg = "ene_c"
  else:
    fnds = "ds_m"
    fndsb = "dsb_m"
    fncfg = "ene_m"
  if not run:
    if iitype == 0: ii = "_aj"
    elif halfwin < 0: ii = "_ix"
    else: ii = "_ii"

    mf = ""
    if halfwin >= 0:
      if mfhalfwin > 0: mf += "_mfii"
      elif mfhalfwin < 0: mf += "_mfav"

    fncfg += ii + mf
    fnds  += ii + mf
    fndsb += ii + mf
  fncfg += ".cfg"
  fnds += ".dat"
  fndsb += ".dat"

  # save configuration file
  open(fncfg, "w").write(strcfg)
  open("ene.cfg", "w").write(strcfg)

  # run
  os.system("../ene")

  # save dat
  shutil.copy("ds.dat", fnds)
  shutil.copy("dsb.dat", fndsb)

''' Main program starts here '''
# build the c program
os.system("make")

# change data dir
os.chdir(datadir)

# canonical
if active:
  run(run = 1)
else:
  shutil.copy2("../data/dsb_c_ii.dat", "dsb.dat")
  shutil.copy2("../data/ds_c_ii.dat", "ds.dat")
run(iitype = 0)
run()
run(mfhalfwin = 70)
run(mfhalfwin = -70)
run(halfwin = -1, mfhalfwin = 30)

# microcanonical ensemble
if active:
  run(canon = 0, run = 1)
else:
  shutil.copy2("../data/dsb_m_ii.dat", "dsb.dat")
  shutil.copy2("../data/ds_m_ii.dat", "ds.dat")
run(canon = 0, halfwin = -1, mfhalfwin = 30)

# call Gnuplot to do the ploting
try:
  import Gnuplot
  g = Gnuplot.Gnuplot(persist = 1)
  g.load(fngp)
  g.reset()
except: # no gnuplot module
  print "cannot find gnuplot module for python, skip plotting"


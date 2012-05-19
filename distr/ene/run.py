#!/usr/bin/env python

''' run the program eh with different conditions, output data '''

import os, sys, re, shutil

active = 0
doextra = 1
datadir = "data"
fngp = "ene.gp"
fncfg = "ene.cfg"

# template for ene_cfg
src_cfg = '''initload = %s
dosimul = %s
nequil = 20000
nstadj = 2000
nsteps = %s 
nstdb = %s # nst of saving to the smaller database dsb
tstat = %s
n = 256
rho = 0.8
rc = 3.0
rs = 2.0
tp = %s
emin = -1500.0
emax = -1100.0
edel = 0.1
iitype = %-4s     # 0: AJ identity; 1: fractional identity
halfwin = %-4s    # 0: single bin; > 0: fixed win; -1: variable window
mfhalfwin = %-4s  # 0: signle bin; > 0: ii; < 0: plain average
# mlimit = 200
sampmin = 400
'''

def run(canon = 1, run = 0, iitype = 1, halfwin = 62, mfhalfwin = 0, save = 1,
    nsteps = 10000000, nstdb = 1000, tp = 1.0):
  # write a config file
  if run:
    initload = 0
    dosimul = 1
  else:
    initload = 1
    dosimul = 0
  strcfg = src_cfg % (initload, dosimul, nsteps, nstdb, canon, tp, iitype, halfwin, mfhalfwin)

  # make file names
  if canon:
    fnds = "ds_c"
    fndsb = "dsb_c"
    fncfg1 = "ene_c"
  else:
    fnds = "ds_m"
    fndsb = "dsb_m"
    fncfg1 = "ene_m"
  if not run:
    if iitype == 0: ii = "_aj"
    elif halfwin < 0: ii = "_ix"
    else: ii = "_ii"

    mf = ""
    if halfwin >= 0:
      if mfhalfwin > 0: mf += "_mfii"
      elif mfhalfwin < 0: mf += "_mfav"

    fncfg1 += ii + mf
    fnds  += ii + mf
    fndsb += ii + mf
  fncfg1 += ".cfg"
  fnds += ".dat"
  fndsb += ".dat"

  # save configuration file
  if save: open(fncfg1, "w").write(strcfg)
  open(fncfg, "w").write(strcfg)

  # run
  os.system("../ene")

  # save dat
  if save: shutil.copy("ds.dat", fnds)
  if save: shutil.copy("dsb.dat", fndsb)

''' Main program starts here '''
# build the c program
os.system("make")

# change data dir
os.chdir(datadir)

# canonical
if active:
  run(run = 1)
else:
  shutil.copy2("../data/dsb_c_ii.dat", "dsb_c.dat")
  shutil.copy2("../data/ds_c_ii.dat", "ds_c.dat")
  shutil.copy2("../data/dsb_c_ii.dat", "dsb.dat")
  shutil.copy2("../data/ds_c_ii.dat", "ds.dat")

run(iitype = 0) # Adib-Jarzynski type 
run()
run(mfhalfwin = 100)
run(mfhalfwin = -100)
run(halfwin = -1, mfhalfwin = 0)

# microcanonical ensemble
if active:
  run(canon = 0, run = 1)
else:
  shutil.copy2("../data/dsb_m_ix.dat", "dsb.dat")
  shutil.copy2("../data/ds_m_ix.dat", "ds.dat")
run(canon = 0, halfwin = -1, mfhalfwin = 0)


if doextra:
  datadirl = "datal"
  datadirh = "datah"
  tlow = 0.8
  thigh = 1.2
  
  # lower temperature
  os.chdir(os.path.join("..", datadirl))
  if active:
    run(run = 1, mfhalfwin = 0, tp = tlow)
  else:
    shutil.copy2("dsb_c.dat", "dsb.dat")
    shutil.copy2("ds_c.dat", "ds.dat")
  run(mfhalfwin = 0, tp = tlow)

  # higher temperature
  os.chdir(os.path.join("..", datadirh))
  if active:
    run(run = 1, mfhalfwin = 0, tp = thigh)
  else:
    shutil.copy2("dsb_c.dat", "dsb.dat")
    shutil.copy2("ds_c.dat", "ds.dat")  
  run(mfhalfwin = 0, tp = thigh)

  # compute error 
  os.chdir("..")
  # compute a version with regular multiple histogram method
  cfg = open(fncfg, "r").readlines()
  cfgbak = cfg[:]
  for i in range(len(cfg)):
    if cfg[i].startswith("halfwin"):
      cfg[i] = "halfwin = 1\n"
      break
  open(fncfg, "w").writelines(cfg)
  os.system("./merge")
  os.chdir(datadir)
  os.rename("dsmerge.dat", "dsmerge0.dat")
  os.rename("dsbmerge.dat", "dsbmerge0.dat")
  os.chdir("..")

  open(fncfg, "w").writelines(cfgbak)
  os.system("./merge")
  os.system("./err > %s/err.dat" % datadir)


  # switch back to data dir  
  os.chdir(datadir)

# call Gnuplot to do the ploting
try:
  import Gnuplot
  g = Gnuplot.Gnuplot(persist = 1)
except: # no gnuplot module
  g = None
  print "cannot find gnuplot module for python, skip plotting %s"  % e

if g:
  os.system("ln -sf ../%s %s" % (fngp, fngp))
  g.load(fngp)
  g.reset()


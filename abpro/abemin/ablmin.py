#!/usr/bin/env python

''' list local minimal energies '''

import os, sys, glob

TMPFILE = "ACCIO.TMP"
verbose = 1

def getsubdirs(pfx):
  ''' list all subdirs under pfx '''
  return sorted( [f for f in os.listdir(pfx) if os.path.isdir(os.path.join(pfx, f))] )

def getop(cmd, err = 0):
  ''' get output from command cmd '''
  if err: err = "2"
  else: err = ""
  os.system(cmd + err + ">" + TMPFILE)
  s = open(TMPFILE).readlines()
  os.remove(TMPFILE)
  return s

def do_datadirs(pfx, verbose):
  ''' list minimal energy from data1 .. '''
  emin = 1e9
  imin = 0
  for i in range(1, 10000):
    dir = "%s/data%s" % (pfx, i)
    fn = os.path.join(dir, "abmin.pos")
    if not os.path.exists(dir): break
    if not os.path.exists(fn): continue
    cmd = "abemin -w -n 10000 -T 1e-12 %s " % fn
    op = getop(cmd, err = 1)
    op = op[0].strip().split()
    ene = float(op[-1])
    if verbose > 1:
      print "%-26s %10.6f" % (fn, ene)
    if ene < emin:
      emin = ene
      imin = i
  if i > 1 and verbose > 0: 
    print "%s/%-6s %10.6f " % (pfx, "data%d"%imin, emin)
  return (imin, emin)

def do_thisdir(dir, verbose):
  ''' do a scratch dir '''
  fn = "%s/abmin.pos" % (dir)
  if not os.path.exists(fn): return (0, 1e10)
  cmd = "abemin -n 100 -T 1e-12 " + fn  + " "
  op = getop(cmd, err = 1)
  op = op[0].strip().split()
  ene = float(op[-1])
  if verbose > 0:
    print "%-26s %10.6f" % (fn, ene)
  return (0, ene)

def do_model(pfx, filter, inscratch, verbose):
  ''' pfx is a model directory like 3dm2, 2dm1, 3dm1 
      filter is 89* or 55* '''
  dirs = sorted(  [f for f in glob.glob(os.path.join(pfx, filter)) if os.path.isdir(f)]  )
  if not len(dirs): return
  #print os.path.abspath(pfx)
  emin = 1e9
  dirmin = "null"
  idmin = 0
  for dir in dirs:
    if inscratch:
      (id, ene) = do_thisdir(dir, verbose - 1)
    else:
      (id, ene) = do_datadirs(dir, verbose - 1)
    if ene < emin:
      emin = ene
      dirmin = dir
      idmin = id
  if dirmin != "null" and verbose > 0:
    if inscratch:
      print "MODEL %s %10.6f\n" % (dirmin, emin)
    else:
      print "MODEL %s/data%s %10.6f\n" % (dirmin, idmin, emin)

def do_models(pfx, inscratch, verbose):
  ''' pfx is a model directory like 3dm2, 2dm1, 3dm1 '''
  do_model(pfx, "34*", inscratch, verbose)
  do_model(pfx, "55*", inscratch, verbose)
  do_model(pfx, "89*", inscratch, verbose)

def do_prj():
  datadirs = getsubdirs(".")
  if "data1" in datadirs:
    do_datadirs(".", verbose)
  else:
    dirs = getsubdirs(".")
    if "3dm2" in dirs:
      for mdl in dirs:
        do_models(mdl, 0, verbose)
    else:
      do_models(".", 0, verbose)

def do_scratch():
  insc = 0
  datadirs = getsubdirs(".")
  dirs = getsubdirs(".")
  if "3dm2" in dirs:
    for mdl in dirs:
      do_models(mdl, 1, verbose)
  else:
    do_models(".", 1, verbose)

def do_args():
  global verbose
  argc = len(sys.argv)
  for i in range(1, argc):
    if sys.argv[i].startswith("-v"):
      verbose = len(sys.argv[i]) # number of v's in -vvv  + 1

''' main function '''
do_args()
path = os.path.abspath(".")
if path.find("project") >= 0:
  do_prj()
else:
  do_scratch()


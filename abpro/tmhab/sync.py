#!/usr/bin/env python

import os, sys

# get hostname
fn = "TMP.OUT"
os.system("hostname > " + fn)
host = open(fn).read().strip()
os.remove(fn)
if host == "tigger": # default to sugar
  server = "cz1@sugar.rice.edu"
else:
  server = "czhang@tigger.rice.edu"

if len(sys.argv) > 1:
  nm =  sys.argv[1]
  if nm == "sugar":
    server = "cz1@sugar.rice.edu"
  elif nm == "tigger":
    server = "czhang@tigger.rice.edu"
  elif nm.startswith("cz5"):
    server = "cz5@sugar.rice.edu"
  else:
    raise Exception

prj = "tmhab"
target = server + ":app/src/" + prj
os.system("rsync -Lvz *.py tmhab.cfg *.[ch] %s" % target)
#os.system("rsync -avz Makefile %s" % target)

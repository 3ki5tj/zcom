#!/usr/bin/env python
''' add 2d-version functions before the 3d-version functions '''


import re, os, sys

sys.path.insert(0, "..")
from python import ccom

ccom.addfd(None, "rc.h", [
  "rc_make", "rc_conj", "rc_norm", "rc_norm2",
  "rc_smul", "rc_sdiv", "rc_add", "rc_sadd", "rc_mul", "rc_div",
  ])


#!/usr/bin/env python
''' add float and double versions of the functions before the real-version functions '''


import re, os, sys

sys.path.insert(0, "..")
from python import ccom

ccom.addfd(None, "rv2.h", [
  "rv2_fprint", "rv2_print",
  "rv2_make", "rv2_makev", "rv2_zero",
  "rv2_copy", "rv2_ncopy", "rv2_swap",
  "rv2_sqr", "rv2_norm", "rv2_dot", "rv2_cross",
  "rv2_neg", "rv2_neg2",
  "rv2_inc", "rv2_dec", "rv2_sinc", "rv2_smul", "rv2_smul2",
  "rv2_normalize", "rv2_makenorm",
  "rv2_diff", "rv2_dist2", "rv2_dist", "rv2_add", "rv2_nadd",
  "rv2_sadd", "rv2_lincomb2",
  "rv2_cosang", "rv2_ang",
  "rv2_rnd", "rv2_rnd0", "rv2_grand",
  "rv2_rnddir0", "rv2_rnddir", "rv2_rnddisk0", "rv2_rnddisk",
  ])

ccom.addfd(None, "rm2.h", [
  "rm2_fprint", "rm2_print",
  "rm2_make", "rm2_makem", "rm2_zero", "rm2_one",
  "rm2_rnduni",
  ])


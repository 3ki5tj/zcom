#!/usr/bin/env python
''' add 2d-version functions before the 3d-version functions '''


import re, os, sys

sys.path.insert(0, "..")
from python import ccom

ccom.addfd(None, "rv3.h", [
  "rv3_fprint", "rv3_print",
  "rv3_make", "rv3_makev", "rv3_zero",
  "rv3_copy", "rv3_ncopy", "rv3_swap",
  "rv3_sqr", "rv3_norm", "rv3_dot", "rv3_cross",
  "rv3_neg", "rv3_neg2",
  "rv3_inc", "rv3_dec", "rv3_sinc", "rv3_smul", "rv3_smul2",
  "rv3_normalize", "rv3_makenorm",
  "rv3_diff", "rv3_dist2", "rv3_dist", "rv3_add", "rv3_nadd",
  "rv3_sadd", "rv3_lincomb2",
  "rv3_cosang", "rv3_ang", "rv3_vdist", "rv3_vpdist",
  "rv3_dih",
  "rv3_rnd", "rv3_rnd0", "rv3_grand", "rv3_grand0",
  "rv3_rnddir0", "rv3_rnddir", "rv3_rndball0", "rv3_rndball",
  ])

ccom.addfd(None, "rm3.h", [
  "rm3_fprint", "rm3_print",
  "rm3_make", "rm3_makem", "rm3_zero", "rm3_one",
  "rm3_copy", "rm3_trans", "rm3_vtv", "rm3_inc", "rm3_sinc",
  "rm3_mul", "rm3_mult", "rm3_tmul", "rm3_mulvec", "rm3_tmulvec",
  "rm3_det", "rm3_inv",
  "rm3_eigval", "rv3_sort3", "rm3_pivot_", "rm3_solvezero",
  "rm3_eigvecs", "rm3_eigsys", "rm3_svd",
  "rm3_mkrot", "rv3_rot",
  "rv3_rmsd",
  "rm3_rnduni",
  ])




#!/usr/bin/env python
'''
add 2d-version functions before the 3d-version functions
currently only on INLINE or static functions
'''


import re, os, sys, copy

sys.path.insert(0, "..")
from python import ccom


def addfd(fnin, fnout, funcls):
  ''' add float and double versions '''

  r2f = {
    "rv3_" : "fv3_",
    "rv2_" : "fv2_",
    "rm3_" : "fm3_",
    "rm2_" : "fm2_",
    "real" : "float",
    }
  r2d = {
    "rv3_" : "dv3_",
    "rv2_" : "dv2_",
    "rm3_" : "dm3_",
    "rm2_" : "dm2_",
    "real" : "double",
    }
  ccom.dupfuncs(fnin, fnout, funcls, [r2f, r2d])


if __name__ == "__main__":
  addfd("rv3.0.h", "rv3.h", [
    "rv3_make", "rv3_zero",
    "rv3_copy", "rv3_ncopy", "rv3_swap",
    "rv3_sqr", "rv3_norm", "rv3_dot", "rv3_cross",
    "rv3_neg", "rv3_neg2",
    "rv3_inc", "rv3_dec", "rv3_sinc", "rv3_smul", "rv3_smul2",
    "rv3_normalize",
    "rv3_makenorm",
    "rv3_diff", "rv3_dist2", "rv3_dist", "rv3_add", "rv3_nadd",
    "rv3_lincomb2",
    "rv3_cosang", "rv3_ang", "rv3_vdist", "rv3_vpdist",
    "rv3_dih",
    "rm3_copy", "rm3_trans", "rm3_vtv", "rm3_inc", "rm3_sinc",
    "rm3_mul", "rm3_mult", "rm3_tmul", "rm3_mulvec", "rm3_tmulvec",
    "rm3_det", "rm3_inv",
    "rm3_eigval", "rv3_sort3", "rm3_pivot_", "rm3_solvezero",
    "rm3_eigvecs", "rm3_eigsys", "rm3_svd",
    "rm3_mkrot", "rv3_rot",
    "rv3_rnd", "rv3_rnd0", "rv3_grand", "rm3_rnduni",
    "rv3_rmsd",
    ])


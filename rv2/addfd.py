#!/usr/bin/env python
'''
add float and double versions of the functions before the real-version functions
currently only on INLINE or static functions
'''


import re, os, sys, copy

sys.path.insert(0, "..")
from python import ccom


def addfd(fnin, fnout, funcls):
  ''' add float and double versions '''

  if fnin == None:
    fnin = ".0".join( os.path.splitext(fnout) )

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
  addfd(None, "rv2.h", [
    "rv2_make", "rv2_zero",
    "rv2_copy", "rv2_ncopy", "rv2_swap",
    "rv2_sqr", "rv2_norm", "rv2_dot", "rv2_cross",
    "rv2_neg", "rv2_neg2",
    "rv2_inc", "rv2_dec", "rv2_sinc", "rv2_smul", "rv2_smul2",
    "rv2_normalize", "rv2_makenorm",
    "rv2_diff", "rv2_dist2", "rv2_dist", "rv2_add", "rv2_nadd",
    "rv2_sadd", "rv2_lincomb2",
    "rv2_cosang", "rv2_ang",
    "rv2_rnd", "rv2_rnd0", "rv2_grand",
    ])



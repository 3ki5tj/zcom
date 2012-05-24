#!/usr/bin/env python
""" create 2d versions from 3d versions """
import os, sys

# add the parent path
sys.path.insert(0, '..')
from python import add2d


if __name__ == "__main__":
  output = "abpro.c"
  s = os.path.splitext(output)
  input = s[0] + ".0" + s[1]
  add2d.add2d(input, output, ["_shake3d", "_rattle3d", "milcshake3d", "milcrattle3d",
      "energy3dm1", "force3dm1", "vv3d"])

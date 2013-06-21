#!/usr/bin/env python
""" create 2d versions from 3d versions """
import os, sys

# add the parent path
sys.path.insert(0, '..')
from python import ccom

ccom.add2d(None, "abpro.h",
    ["_shake3d", "_rattle3d", "milcshake3d", "milcrattle3d",
     "energy3dm1", "force3dm1", "vv3d"])


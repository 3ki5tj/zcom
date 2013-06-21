#!/usr/bin/env python
""" create 2d versions from 3d versions """
import os, sys

# add the parent path
sys.path.insert(0, '..')
from python import ccom

ccom.add2d(None, "lj.h", [
  "lj_pbcdist2_3d",
  "lj_energysw3d", "lj_energysq3d", "lj_energylj3d", "lj_energyx3d", "lj_energy3d",
  "lj_forcesw3d", "lj_forcelj3d", "lj_force3d",
  ])

ccom.add2d(None, "ljmc.h", [
  "lj_randmv3d",
  "lj_pairsq3d", "lj_depotsq3d", "lj_commitsq3d", "lj_metrosq3d",
  "lj_pairsw3d", "lj_depotsw3d", "lj_commitsw3d", "lj_metrosw3d",
  "lj_pairlj3d", "lj_depotlj3d", "lj_commitlj3d", "lj_metrolj3d",
  "lj_pair3d", "lj_depot3d", "lj_commit3d", "lj_metro3d",
  "lj_dupertl3d", "lj_dupertg3d", "lj_duinsert3d",
  ])

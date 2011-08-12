#!/usr/bin/env python

"""
create 2d versions from 3d versions
"""

import re, os, copy


def add2d(lines, funcname):
  ''' add 2D version of funcname '''

  # 1. search the beginning
  for i0 in range(len(lines)):
    s = lines[i0]
    if s.startswith("static") and s.find(funcname+"(") >= 0:
      break
  else:
    print "cannot find function start of %s" % funcname
    raise Exception

  # 2. search the end
  for i1 in range(i0+1, len(lines)):
    if lines[i1].rstrip() == '}':
      break
  else:
    print "cannot find function end"

  func3 = lines[i0: i1+1]
  func2 = ["\n"]*len(func3)
  for i in range(len(func2)):
    s = func3[i]
    s = re.sub("rv3_", "rv2_", s)
    s = re.sub(r"\[3\]", "[2]", s)
    func2[i] = s
  func2name = re.sub("3d", "2d", funcname)
  func2[0] = re.sub(funcname, func2name, func2[0])

  lines = lines[:i0] + func2 + ["\n"] + lines[i0:]
  return lines


if __name__ == "__main__":
  output = "abpro.c"
  s = os.path.splitext(output)
  input = s[0] + ".0" + s[1]
  s = open(input, 'r').readlines()
  for f in ["_shake3d", "_rattle3d", "_milcshake3d", "energy3dm1", "force3dm1"]:
    s = add2d(s, f)
  s = ["/* Don't edit this file, edit %s instead. */\n" % input,] + s
  open(output, 'w').writelines(s)


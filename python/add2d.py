#!/usr/bin/env python
'''
add 2d-version functions before the 3d-version functions
currently only on INLINE or static functions
'''


import re, os, copy

def addfunc2d(lines, funcname):
  ''' add 2D version of a 3D function '''

  ismacro = 0

  # 1. search the beginning of the function
  for i0 in range(len(lines)):
    s = lines[i0]
    # we only operate on static and inline functions
    if ( (s.startswith("static") or s.startswith("INLINE") )
        and s.find(funcname + "(") >= 0):
      break
    if ( s.startswith("#define")
        and s.find(funcname + "(") >= 0):
      ismacro = 1
      break
  else:
    print "cannot find function start of %s" % funcname
    raise Exception

  # 2. search the end of the function
  for i1 in range(i0, len(lines)):
    if ismacro: # macro
      if lines[i1].strip() == "" or lines[i1].strip()[-1] != "\\":
        break
    else:  #  function
      if lines[i1].rstrip() == '}':
        break
    s = lines[i1].strip()
    if s.startswith("{") and s.endswith("}") and lines[i1 + 1].strip() == "":
      break
  else:
    print "cannot find function end"
    raise Exception

  # 3. move upwards until the comment starts
  while i0 > 0:
    if lines[i0 - 1].strip() == "": break
    if lines[i0].strip().startswith("/*"): break
    i0 -= 1

  func3 = lines[i0 : i1 + 1]
  func2 = ["\n"] * len(func3)
  for i in range(len(func2)):
    s = func3[i]
    s = re.sub("rv3_", "rv2_", s)
    s = re.sub("rm3_", "rm2_", s)
    s = re.sub("dv3_", "dv2_", s)
    s = re.sub("dm3_", "dm2_", s)
    s = re.sub("fv3_", "fv2_", s)
    s = re.sub("fm3_", "fm2_", s)
    s = re.sub(r"\[3\]", "[2]", s)
    s = re.sub(r"3d\(", "2d(", s)
    s = re.sub(r"d < 3", "d < 2", s)
    s = re.sub(r"3(\s*)\*", r"2\1*", s)
    s = re.sub(r"\*3", r"*2", s)
    s = re.sub(r"\* 3", r"* 2", s)
    s = re.sub(r"3D", "2D", s)
    func2[i] = s
  func2name = re.sub("3d", "2d", funcname)
  func2[0] = re.sub(funcname, func2name, func2[0])

  lines = lines[:i0] + func2 + ["\n"] + lines[i0:]
  return lines

def add2d(fnin, fnout, funcls):
  ''' add 2D version '''
  s = open(fnin, "r").readlines()
  for f in funcls:
    s = addfunc2d(s, f)
  s = ["/* Don't edit this file, edit %s instead. */\n"
       % fnin,] + s
  if fnout:
    open(fnout, 'w').writelines(s)
  else:
    print ''.join(s)

if __name__ == "__main__":
  add2d("add2dtest.c", None, ["func3d", ])


#!/usr/bin/env python

''' common routines for manipulating C source code '''

import re



def getfunc(lines, funcname):
  ''' return the begin and ending of a function '''

  ismacro = 0

  # 1. search the beginning of the function
  for i0 in range(len(lines)):
    s = lines[i0]
    # we only operate on static and inline functions
    if ( (s.startswith("static") or s.startswith("INLINE") )
        and funcname + "(" in s:
      break
    if ( s.startswith("#define")
        and funcname + "(" in s:
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

  return i0, i1 + 1



def dictrepl(lines, d):
  ''' make replacements according to the dictionary `d'
      keys and values in the dictionary are regular expressions '''

  n = len(lines)
  dest = [""] * n
  for i in range(n):
    s = lines[i]
    for k in d:
      s = re.sub(k, d[k], s)
    dest[i] = s
  return dest



def dupfunc(lines, funcname, d):
  ''' duplicate function `funcname' with modifications
      given by the dictionary `d' '''

  i0, i1 = getfunc(lines, funcname)
  func2 = dictrepl(lines[i0 : i1], d)
  return lines[:i0] + func2 + ["\n"] + lines[i0:]



def dupfuncs(fnin, fnout, funcls, dls):
  ''' duplicate all functions in funcls with modifications
      given by the dictionary list `dls' '''

  # make a dictionary list, if the input is a single dictionary
  if type(dls) == dict:  dls = [ dls, ]

  s = open(fnin, "r").readlines()
  for f in funcls:
    for d in dls:
      s = dupfunc(s, f, d)
  s = ["/* Don't edit this file, edit %s instead. */\n"
       % fnin,] + s
  if fnout:
    open(fnout, 'w').writelines(s)
  else:
    print ''.join(s)



def add2d(fnin, fnout, funcls):
  ''' add 2D versions '''

  dupfuncs(fnin, fnout, funcls, {
      "rv3_" : "rv2_",
      "rm3_" : "rm2_",
      "fv3_" : "fv2_",
      "fm3_" : "fm2_",
      "dv3_" : "dv2_",
      "dm3_" : "dm2_",
      r"\[3\]" : "[2]",
      r"3d\(" : "2d(",
      r"d < 3" : "d < 2",
      r"3(\s*)\*" : r"2\1*",
      r"\*3" : r"*2",
      r"\* 3": r"* 2",
      r"3D" : "2D",
      r"3d" : "2d",
    })


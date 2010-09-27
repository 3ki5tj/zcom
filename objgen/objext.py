#!/usr/bin/env python
'''
additional components 
'''

def notalways(cond):
  ''' test if a condition is missing or always true '''
  return cond not in (None, 1, "1", "TRUE")

class Fold:
  ''' a sub-object embedded inside an object '''
  def __init__(f, fprefix, name):
    # generate a name
    f.name = name
    # generate a function prefix
    if not name.endswith("_"):
      namep = name + ("_" if len(name) else "")
    f.fprefix = fprefix + namep
    f.items = []

  def additem(f, item):
    f.items += [item]

  def __len__(f):
    return len(f.items)

  def add_fold_tests(f, cow, funcnm, ret = "0", validate = 1):
    ''' add prerequisite/validation tests '''
    if notalways(f.prereq):
      cow.addln("if ( !(%s) ) return %s;", f.prereq, ret)
    if notalways(f.valid) and validate:
      cow.validate(f.valid, funcnm)


def type2fmt_s(tp):
  ''' format for scanf '''
  fmt = ""
  if   tp == "int":      
    fmt = "%d"
  elif tp in ("long", "long int"):
    fmt = "%ld"
  elif tp in ("long long", "long long int"):
    fmt = "%lld"
  elif tp in ("unsigned", "unsigned int"): 
    fmt = "%u"
  elif tp == "unsigned long": 
    fmt = "%ul"
  elif tp == "float":
    fmt = "%f"
  elif tp == "double":
    fmt = "%lf"
  elif tp == "char *":
    fmt = "%s"
  else:
    print "no scanf format string for type [%s]" % (tp)
    raise Exception
  return fmt

def type2fmt_p(tp, var):
  ''' format for printf '''
  fmt = ""
  if   tp in ("int", "char"):      
    fmt = "%4d"
  elif tp in ("long", "long int"):
    fmt = "%4ld"
  elif tp in ("long long", "long long int"):
    fmt = "%8lld"
  elif tp in ("unsigned", "unsigned int"): 
    fmt = "0x%8X"
  elif tp == "unsigned long": 
    fmt = "0x%8lX"
  elif tp in ("float", "double"):
    fmt = "%f"
  elif tp == "char *":
    fmt = "%s"
  else:
    #print "no printf format string for type [%s]%s" % (tp,
    #    " item: "+var if var else "")
    raise TypeError
  return fmt

def checkcycl_i(deps, i0, n, checked):
  ''' check cyclic dependencies starting from i0 '''
  checked[i0] = 1
  st = [0] * (n+1)
  jt = [0] * (n+1)
  top = 0
  st[top] = i0
  jt[top] = 0
  while top >= 0:
    i = st[top]
    j = jt[top]
    dl = deps[i] # list of i depends on
    if j >= len(dl):
      top -= 1
      if top >= 0:
        jt[top] += 1
    else:
      top += 1
      x = dl[j]
      if x in st[:top]:
        print "cyclic dependence is detected from %s, stack %s" % (x, st[:top])
        return -1
      st[top] = x
      jt[top] = 0
      checked[x] = 1 # circular dependency from x will be checked as well
    # print "i0: %d, top: %s, stack: %s, j: %s" % (i0, top, st[:top+1], jt[:top+1]); raw_input()
  return 0

def checkcycl(deps):
  ''' check circular dependencies '''
  n = len(deps)
  checked = [0] * n;
  for i in range(n):
    if checked[i]: continue
    if 0 != checkcycl_i(deps, i, n, checked):
      return -1
  return 0

def sortidx(deps, verbose = 0):
  ''' 
  return a sorted index array, depending on a dependency array
  using a bubble sort
  `deps' is not changed 
  '''
  n = len(deps)
  d = range(n)  # dictionary
  i = 0
  if checkcycl(deps) != 0:
    print "dependency is cyclic!"
    return None
  while i < n:
    if verbose: print "i: %d -> %d  deps: %s\n%s" % (i, d[i], deps[d[i]], d)
    if len(deps[d[i]]) == 0: 
      i += 1
      continue
    if verbose: print "i %d, actual deps: %s" % (i, map(d.index, deps[d[i]]) ); raw_input()
    mi = max(map(d.index, deps[d[i]]))
    if mi <= i:
      i += 1
      continue
    # swap i with mi
    d[i], d[mi] = d[mi], d[i]
  return d

def testrank():
  n = 8
  deps = [[]] * n
  deps[0] = [6]
  deps[1] = [2, 3]
  deps[2] = [4]
  deps[3] = [4]
  deps[4] = []
  deps[5] = [7]
  #checkcycl(deps)
  print deps
  print sortidx(deps, 0)

if __name__ == "__main__":
  testrank()



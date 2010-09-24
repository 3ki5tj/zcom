#!/usr/bin/env python
'''
additional components 
'''

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

def dump_block(block, cw, tab = 4, offset = 2):
  ''' 
  dump a block of lines that are broken to fields 
  cw is a CodeWriter
  '''
  if len(block) == 0: return
  nfields = len(block[0])
  wid = [8] * nfields
  #print "block of %s\n%s" % (len(block), block)
  #raw_input()

  # align and format all items in the bufferred block
  # and append the result to string s
  for j in range(nfields):
    w = max(len(item[j]) for item in block) + 1 # +1 for the following blank
    wid[j] = ((w + offset + tab - 1) // tab)*tab - offset

  for item in block:
    s = ""
    for j in range(nfields-1):
      s += "%-*s" % (wid[j], item[j])
    s += item[nfields - 1] 
    cw.addln(s.rstrip())

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



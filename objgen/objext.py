#!/usr/bin/env python
'''
additional components 
'''
import re

MASTERID = "0"       # master rank in MPI

def isint(s):
  ''' judge if a string is integer '''
  if type(s) == int:
    return 1
  try:
    i = int(s)
  except ValueError:
    return 0
  return 1

def notalways(cond):
  ''' test if a condition is missing or always true '''
  return cond not in (None, 1, "1", "TRUE")

def disabled(cond):
  ''' 
  test if a default command is explicitly disabled
  cond is assumed to be enabled if it is None
  '''
  return cond in (0, "0", "FALSE", "off")

def findvar(var, s):
  pattern = "(?<![a-z_])%s(?![a-z_])" % var
  m = re.search(pattern, s)
  return (m.start(0) + 1) if m else 0

def escape(s, printf=1):
  if type(s) != str:
    print "s is not string! s = %s" % s
    raise Exception
  s = re.sub(r'"', r'\"', s)
  s = re.sub(r"'", r"\'", s)
  if printf:
    s = re.sub("\%", "%%", s)
  return s

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
    fmt = "0x%X"
  elif tp == "unsigned long": 
    fmt = "0x%lX"
  elif tp in ("float", "double", "real"):
    fmt = "%g"
  elif tp == "char *":
    fmt = "%s"
  elif tp.endswith("pointer"):
    fmt = "%p"
  else:
    #print "no printf format string for type [%s]%s" % (tp,
    #    " item: "+var if var else "")
    raise TypeError
  return fmt

def type2zero(tp):
  zval = None
  if tp in ("int", "unsigned", "unsigned int", 
      "long", "unsigned long"):
    zval = "0"
  elif tp in ("float", "real"):
    zval = "0.0f"
  elif tp in ("double", "long double"):
    zval = "0.0"
  elif tp in ("pointer", "function pointer", "char *", 
      "pointer to object"):
    zval = "NULL"
  elif tp == "MPI_Comm":
    zval = "MPI_COMM_NULL"
  return zval

def mpitype(tp):
  if tp == "double": return "MPI_DOUBLE"
  elif tp == "float": return "MPI_FLOAT"
  elif tp == "int": return "MPI_INT"
  elif tp == "long": return "MPI_LONG"
  elif tp == "unsigned char": return "MPI_UNSIGNED_CHAR"
  elif tp == "byte": return "MPI_BYTE"
  else: raise Exception

def remove_idle_pp(lines):
  ''' remove empty preprocessor blocks 
    #if
    #else
    #endif
  similar to merge_pp_if_blocks
  '''
  i = 1
  # remove empty pp
  while i < len(lines):
    if lines[i].startswith("#endif") and i > 0:
      lnprev = lines[i-1]
      if (lnprev.startswith("#else") or 
          lnprev.startswith("#elif")):
        # remove an empty #else/#elif branch
        lines = lines[:i-1] + lines[i:]
        continue
      elif lnprev.startswith("#if"):
        lines = lines[:i-1] + lines[i+1:]
        continue
    i += 1

  # merge neighoring preprocessing block
  # if their conditions are the same
  cond = None
  i = 1
  while i < len(lines):
    line = lines[i].strip()
    if line.startswith("#if "):
      if not cond: cond = line[4:].strip()
    elif line.startswith("#endif"):
      if (not cond or i >= len(lines) - 2
          or not lines[i+1].startswith("#if ")):
        cond = None
        i += 1
        continue
      ncond = lines[i+1].strip()[4:].strip()
      if cond == ncond:
        # print "removing\n%s\n%s" % (lines[i], lines[i+1]); raw_input()
        lines = lines[:i] + lines[i+2:]
      else:
        cond = None
    elif line.startswith("#el"):
      cond = None
    i += 1
  return lines

def next_block(lines, i, cond):
  ''' return the line index if after several blank lines
  a line exactly equal to cond occurs, otherwise return -1 '''
  n = len(lines)
  for j in range(i+1, n):
    if lines[j].strip() != "":
      break
  else:
    return -1
  if j >= n-1: return -1
  if lines[j] == cond:
    return j
  else: return -1

def merge_if_blocks(lines):
  '''
  merge neighboring if-blocks with the same condition
      if (abc) {
        ...
  X   }
  X   if (abc) {
        ...
      }
  '''
  return merge_blocks(lines, r"\s*if\s*\(.*\)\s*{$", "if", r"\s*(\})\s*else", "}")

def merge_pp_if_blocks(lines):
  return merge_blocks(lines, r"\s*#if.*$", "#if", "\s*(#el)", "#endif")

def merge_blocks(lines, pattern, sbegin, selse, send):
  cond = None
  i = 1
  while i < len(lines):
    line = lines[i]

    # if
    m = re.match(pattern, line.strip())
    if m:
      #print "i:%4d, %s" % (i, line); raw_input()
      if not cond: cond = line
      i += 1
      continue

    # else
    m = re.match(selse, line)
    if m:
      if cond and m.start(1) <= cond.find(sbegin):
        cond = None
      i += 1
      continue

    # finish an if block
    if line.strip() == send and (cond and 
        cond.find(sbegin) >= line.find(send) ):
      if cond != None and cond.find(sbegin) == line.find(send):
        ip = next_block(lines, i, cond)
        #if sbegin == "if":
        #  print "match end i=%d ip=%d\n%s\n%s\n" % (i, ip, cond, line); raw_input()
        if ip >= 0:
          #print "removing i: %d - %d\n%s\n%s\n...%s\n" % (i, ip, lines[i], lines[i+1], lines[i+2]); raw_input()
          lines = lines[:i] + lines[i+1:ip] + lines[ip+1:]
          continue
      cond = None  # terminate condition
    i += 1
  return lines

def use_ifdef(lines):
  pat = r"\s*#if\s+(\!?)defined\((\w+)\)$"
  n = len(lines)
  for i in range(n):
    line = lines[i]
    #print "(%s)"%line; raw_input()
    m = re.match(pat, line)
    if not m: continue
    pfx = "#ifndef" if len(m.group(1)) else "#ifdef"
    lines[i] = pfx + " " + m.group(2)
    #print "[%s] -> [%s]" % (line, lines[i]); raw_input()
  return lines

def trimcode(s):
  ''' various code simplification '''
  endl = '\n' if s.endswith('\n') else ''
  lines = s.splitlines()
  lines = merge_pp_if_blocks(lines)
  lines = remove_idle_pp(lines)
  lines = use_ifdef(lines)
  lines = merge_if_blocks(lines)
  return '\n'.join(lines) + endl

# item sorting

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



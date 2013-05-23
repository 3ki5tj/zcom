#!/usr/bin/env python

fn = "spgr.h"

import os, sys, re

src = open(fn).readlines()
n = len(src)

output = ""


# read types
types = ""
for i in range(n):
  if src[i].find("spacegroup_types") >= 0:
    break
else:
  raise Exception

for j in range(i + 1, i + 532):
  ln = src[j].strip()
  ln = ln.replace("{", "")
  ln = ln.replace("}", "")
  while 1:
    m = re.search(r'"\s*(.*?)\s*"', ln)
    if not m: break
    ln = ln[:m.start(0)] + '$' + m.group(1).strip() + '$' + ln[m.end(0):]
  ln = ln.replace("$", '"')
  k = ln.find("/*")
  ln = ln[:k].strip()
  if ln.endswith(","): ln = ln[:-1]
  arr = ln.split(",")
  arr = arr[:-1]
  ln = "[" + ",".join( arr ) + "],\n"
  types += ln
types = "spacegroup_types = [\n" + types + "];\n\n\n"

output += types

# do symmetry_operations
symops = ""
for i in range(n):
  if src[i].find("symmetry_operations") >= 0:
    break
else:
  raise Exception

for j in range(i + 1, i + 7390):
  ln = src[j].strip()
  k = ln.find(",")
  ln = ln[:k].strip() + ",\n"
  symops += ln
symops = "symmetry_operations = [\n" + symops + "];\n\n\n"

output += symops

# read indices
indices = ""
for i in range(n):
  if src[i].find("symmetry_operation_index[][2]") >= 0:
    break
else:
  raise Exception

for j in range(i + 1, i + 532):
  ln = src[j].strip()
  ln = ln.replace("{", "")
  ln = ln.replace("}", "")
  k = ln.find("/*")
  ln = ln[:k].strip()
  if ln.endswith(","): ln = ln[:-1]
  arr = ln.split(",")
  arr = [s.strip() for s in arr]
  ln = "[" + ",".join( arr ) + "],\n"
  indices += ln
indices = "symmetry_operation_index = [\n" + indices + "];\n\n\n"

#print indices
output += indices


head = """#!/usr/bin/env python

import os, sys, re

"""


tail = """

def getop(id):
  w = symmetry_operations[id]
  r = w % (3**9)
  d = 3**8
  rot = [[0,0,0],[0,0,0],[0,0,0]];
  for i in range(3):
    for j in range(3):
      rot[i][j] = (r % (d * 3)) / d - 1
      d /= 3

  t = w / (3**9)
  d = 144
  trans = [0,0,0]
  for i in range(3):
    trans[i] = 1.*((t % (d*12))/d)/12;
    d /= 12;
  return rot, trans

def getmat(hall):
  lsm = []
  lsv = []
  cnt = symmetry_operation_index[hall][0]
  ist = symmetry_operation_index[hall][1]
  for i in range(cnt):
    rot, trans = getop(ist + i)
    lsm += [rot]
    lsv += [trans]
  return lsm, lsv

def iscsym(hall):
  cnt = symmetry_operation_index[hall][0]
  ist = symmetry_operation_index[hall][1]
  for i in range(cnt):
    if symmetry_operations[ist + i] == 3198: # inversion
      return 1
  return 0

def getmatstr(s, offset = 0):
  ''' get spacegroup from string '''
  SPMAX = 530

  for i in range(1, SPMAX+1):
    buf1 = spacegroup_types[i][3].replace(" ", "")
    buf2 = spacegroup_types[i][4].replace(" ", "")
    buf3 = spacegroup_types[i][5].replace(" ", "")
    if (s == buf1 or s == buf2 or s == buf3):
      break
  else:
    print "no match for %s" % s
    return None, None, 0
  i += offset
  lsm, lsv = getmat(i)
  return lsm, lsv, i


if __name__ == "__main__":
  spgr = "P4_122"
  if len(sys.argv) > 1:
    spgr = sys.argv[1]
  ms, vs, n = getmatstr(spgr)
  for i in range(len(ms)):
    print ms[i], vs[i]
  print spacegroup_types[n][0]

  cnt = len(spacegroup_types)
  ls = []
  for i in range(1, cnt):
    tp = spacegroup_types[i]
    csym = iscsym(i)
    if not csym and not tp[5] in ls:
      ls = ls + [ tp[5] ]
  print ls
"""

output = head + output + tail

open("spgr_draft.py", "w").write(output)


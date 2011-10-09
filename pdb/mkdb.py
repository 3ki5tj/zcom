#!/usr/bin/env python

import copy

''' generate indices for proteins '''

aa_atoms_str = '''
GLY HA1 HA2  # HA3
ALA CB HB1 HB2 HB3
VAL CB HB        CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23
LEU CB HB1 HB2   CG  HG   CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23
ILE CB HB        CG2 HG21 HG22 HG23 CG1 HG11 HG12 CD HD1 HD2 HD3
PRO CD HD1 HD2   CG  HG1  HG2  CB HB1 HB2 # HG3
SER CB HB1 HB2 OG HG
THR CB HB  CG2 HG21 HG22 HG23 OG1 HG1
CYS CB HB1 HB2 SG HG
MET CB HB1 HB2 CG HG1 HG2 SD CE HE1 HE2 HE3 # HG3
ASN CB HB1 HB2 CG OD1 ND2 HD21 HD22
GLN CB HB1 HB2 CG HG1 HG2 CD OE1 NE2 HE21 HE22 # HG3
ASP CB HB1 HB2 CG OD1 OD2
GLU CB HB1 HB2 CG HG1 HG2 CD OE1 OE2 # HG3
LYS CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 CE HE1 HE2 NZ HZ1 HZ2 HZ3 # HG3
ARG CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 # HG3
HIS CB HB1 HB2 CG ND1 CE1 HE1 NE2 HE2 CD2 HD2
PHE CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2
TYR CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ OH HH CE2 HE2 CD2 HD2
TRP CB HB1 HB2 CG CD1 HD1 NE1 HE1 CE2 CZ2 HZ2 CH2 HH2 CZ3 HZ3 CE3 HE3 CD2
'''.strip().upper()

ss = aa_atoms_str.split("\n")
if len(ss) != 20: raise Exception
atnm = [""]*20
nmdb = []
for i in range(20):
  # including alternative names
  k = ss[i].find("#")
  line = ss[i] if k < 0 else ss[i][:k] + ss[i][k+1:]
  
  line = line.split()
  aanm = line[0]
  atnm[i] = line[1:]
  atnm[i] += "CA N H2 H3 C OC1 OC2".strip().split()
  if aanm != "GLY": atnm[i] += [ "HA" ]
  if aanm != "PRO": atnm[i] += [ "H1" ]
  for s in atnm[i]:
    if s not in nmdb:
      nmdb += [ s ]

#print len(nmdb), "items: ", ' '.join(nmdb)

users = [ "" ] * len(nmdb)
for k in range(len(nmdb)):
  nm = nmdb[k]
  users[k] = [nm, ] # first element is the atom name
  for i in range(20):
    if nm in atnm[i]:
      users[k] += [ i ]
  #print "%s %s" % (nm, users[k])

# sort the atom name array by # of users
users = sorted(users, key = lambda x: -len(x))
#print "%d users" % (len(users))
#users1 = copy.deepcopy(users)

def get_bitmap(users):
  bitmap = {}
  bit = 0
  while bit < len(users):
    ent = users[bit]
    at = ent[0]
    usr = ent[1:]
    set1 = set(usr)

    # search for largest non-overlapping set
    jmax = -1
    szmax = 0
    for j in range(i+1, len(users)):
      set2 = set(users[j][1:])
      intsect = set1 & set2
      nol = len(intsect)
      if nol > 0: continue # overlap
      sz = len(set2)
      if sz > szmax:
        szmax = sz
        jmax = j
    if jmax < 0:
      #print "no non-overlapping atom for %s, covers %d" % (at, len(set1))
      bitmap[at] = bit
      bit += 1
      continue
    else:
      at2 = users[jmax][0]
      #print "largest non-overlapping atom for %s is %s" % (at, at2)
      bitmap[at] = bit
      bitmap[at2] = bit
      # don't increment bit, because, we're still working on growing the set
      users[bit] += users[jmax][1:] # add users from jmax
      users = users[:jmax] + users[jmax+1:] # remove jmax
      #print "remove %s at %d, new users %d: %s, len = %d" % (at2, jmax, bit, users[bit], len(users))
    #raw_input()
  # create backbone aliases
  bitmap['H'] = bitmap['H1']
  bitmap['O'] = bitmap['OC1']
  bitmap['OXT'] = bitmap['OC2']
  return sorted(bitmap.iteritems(), key = lambda (k,v): (v,k))

# sort by bit index
bitmap = get_bitmap(users)

def declare_bitmap(bitmap, pfx):
  ''' C declaration '''
  for (k, v) in bitmap:
    print "#define %-9s %s" % (pfx+k, v)

prefix = "AA_"

declare_bitmap(bitmap, prefix)

def write_code(bitmap, pfx):
  str = 'static struct { unsigned bit; const char *nm; } aaatomtypes[] = {\n  '
  for i in range(len(bitmap)):
    (nm, bit) = bitmap[i]
    str += '{%-8s %-6s}, ' % (pfx+nm+",", '"'+nm+'"')
    if  (i+1) % 5 == 0:
      str += '\n  '
  str = str.strip() + '\n  {-1, NULL}};\n';
  print str
write_code(bitmap, prefix)

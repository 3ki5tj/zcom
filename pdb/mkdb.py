#!/usr/bin/env python

import copy

''' generate indices for proteins '''

dbg = 2

aa_atoms_str = '''
GLY HA1 HA2  # HA3:HA1
ALA CB HB1 HB2 HB3
VAL CB HB        CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23
LEU CB HB1 HB2   CG  HG   CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23 # HB3:HB1
ILE CB HB        CG2 HG21 HG22 HG23 CG1 HG11 HG12 CD HD1 HD2 HD3 # HG13:HG11 CD1:CD HD11:HD1 HD12:HD2 HD13:HD3
PRO CD HD1 HD2   CG  HG1  HG2  CB HB1 HB2 # HB3:HB1 HG3:HG1 HD3:HD1
SER CB HB1 HB2 OG HG # HB3:HB1
THR CB HB  CG2 HG21 HG22 HG23 OG1 HG1
CYS CB HB1 HB2 SG HG # HB3:HB1
MET CB HB1 HB2 CG HG1 HG2 SD CE HE1 HE2 HE3 # HB3:HB1 HG3:HG1
ASN CB HB1 HB2 CG OD1 ND2 HD21 HD22 # HB3:HB1
GLN CB HB1 HB2 CG HG1 HG2 CD OE1 NE2 HE21 HE22 # HB3:HB1 HG3:HG1
ASP CB HB1 HB2 CG OD1 OD2 # HB3:HB1
GLU CB HB1 HB2 CG HG1 HG2 CD OE1 OE2 # HB3:HB1 HG3:HG1
LYS CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 CE HE1 HE2 NZ HZ1 HZ2 HZ3 # HB3:HB1 HG3:HG1 HD3:HD1 HE3:HE1
ARG CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 # HB3:HB1 HG3:HG1 HD3:HD1
HIS CB HB1 HB2 CG ND1 HD1 CE1 HE1 NE2 HE2 CD2 HD2 # HB3:HB1
PHE CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2 # HB3:HB1
TYR CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ OH HH CE2 HE2 CD2 HD2 # HB3:HB1
TRP CB HB1 HB2 CG CD1 HD1 NE1 HE1 CE2 CZ2 HZ2 CH2 HH2 CZ3 HZ3 CE3 HE3 CD2 # HB3:HB1
'''.strip().upper()

def main():
  ss = aa_atoms_str.split("\n")
  if len(ss) != 20: raise Exception
  resnm, resat, resub = parse_atom_names(ss)
  declare_database(resnm, resat, resub)

def parse_atom_names(ss):
  ''' returns
  resnm: array of residues
  resat: array of atom names in each residue
  '''
  resnm = [""]*20
  resat = [""]*20  # list of atom names in each residue
  resub = [""]*20
  for i in range(20):
    # including alternative names
    k = ss[i].find("#")
    line = ss[i] if k < 0 else ss[i][:k]
    subs = "" if k < 0 else ss[i][k+1:]
    line = line.split()
    resnm[i] = aanm = line[0]
    resat[i] = "CA N C O".split() + line[1:]
    resub[i] = subs.strip().split()
    if aanm != "GLY": resat[i] += [ "HA" ]
    if aanm != "PRO": resat[i] += [ "H" ]
  return resnm, resat, resub

def declare_database(resnm, resat, resub):
  maxatm = max([len(ent) for ent in resat]) + 1
  maxsub = max([len(ent) for ent in resub])*2 + 1
  strtag = "tag_pdb_aadb"
  ss = '''
struct %s {
  const char *resnm;      /* residue name */
  const char *atom[%s];   /* atoms */
  const char *sub[%s];    /* substitutions */
  unsigned long hvflags;  /* backbone and heavy atom flags */
} pdb_aadb[20] = {'''.lstrip() % (strtag, maxatm, maxsub) + "\n"

  for i in range(20):
    s = '{"%s", ' % resnm[i]
    # list atoms
    s += '{'
    for at in resat[i]: s += '"%s", ' % at
    s += 'NULL}, '

    # list substitutions
    s += '{'
    for sub in resub[i]:
      alt = sub.split(":")
      s += '"%s", "%s", ' % (alt[0], alt[1])
    s += 'NULL}, '

    # heavy
    hvflags = 0
    for k in range(len(resat[i])):
      at = resat[i][k]
      if at[0] != 'H':
        hvflags |= 1 << k
    s += "0x%xul" % hvflags
    s += "},\n"
    ss += s
  ss = ss.rstrip()[:-1] + "};\n"
  #print ss

  fn = "pdb.h"
  src = open(fn).readlines()
  for i in range(len(src)):
    if src[i].strip().startswith("struct %s" % strtag):
      break
  else:
    print "Error: cannot find insertion point!\n%s" % ss
    raise Exception
  for j in range(i, len(src)):
    if src[j].strip() == "": break
  src = src[:i] + [ ss ] + src[j:]
  open(fn, "w").writelines(src)

# TODO: handle H1, H2, H3, OXT, OC1 OC2


if __name__ == "__main__":
  main()


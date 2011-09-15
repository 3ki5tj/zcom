#ifndef PDB_H__
#define PDB_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct {
  real x[3];
  int aid; /* atom index */
  int rid; /* residue index */
  int insert;
  char atnm[8];
  char resnm[8];
} pdbatom_t; /* a single atom entry */

typedef struct {
  int n; /* number of lines */
  int nalloc;
  int nres;
  pdbatom_t *at;
  const char *file;
} pdbmodel_t; /* raw data in a pdb model */

#define HAS_N   1
#define HAS_CA  2
#define HAS_C   4
#define HAS_O   8
#define HAS_BB  (HAS_N|HAS_CA|HAS_C|HAS_O)
#define HAS_CB 16

typedef struct {
  real xn[3], xca[3], xc[3], xo[3], xcb[3];
  int aa;
  int broken;
  unsigned flags;
} pdbaar_t; /* amino-acid residues */

typedef struct {
  pdbaar_t *res;
  int nres;
  const char *file;
} pdbaac_t; /* amino-acid chain */

/* generic pdb model */
pdbmodel_t *pdbm_read(const char *fname, int verbose);
int pdbm_write(pdbmodel_t *m, const char *fn);
#define pdbm_free(m) { free(m->at); free(m); }

/* protein pdb */
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose);
#define pdbaac_free(c) { free(c->res); free(c); }

static const char pdb_aanames_[21][4] = {"---",
  "GLY", "ALA", "VAL", "LEU", "ILE", "PRO",
  "THR", "SER", "CYS", "MET", "PHE", "TYR", "TRP",
  "GLU", "GLN", "ASP", "ASN", "ARG", "LYS", "HIS"};

ZCINLINE int pdbaaidx(const char *res)
{
  int i;
  for (i = 1; i <= 20; i++)
    if (strcmp(res, pdb_aanames_[i]) == 0) 
      return i;
  return -1;
}

ZCINLINE const char *pdbaaname(int i)
{
  if (i <= 0 || i > 20) {
    fprintf(stderr, "invalid pdb id %d\n", i);
    exit(1);
  }
  return pdb_aanames_[i];
}


#endif


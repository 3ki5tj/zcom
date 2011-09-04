#ifndef PDB_H__
#define PDB_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct {
  real x[3];
  int aid;
  int rid;
  int insert;
  char atnm[8];
  char resnm[8];
} pdbatom_t;

typedef struct {
  int n; /* number of lines */
  int nalloc;
  int nres;
  pdbatom_t *at;
  const char *file;
} pdbmodel_t;

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
} pdbaares_t;

typedef struct {
  pdbaares_t *res;
  int nres;
  const char *file;
} pdbaabb_t;

pdbmodel_t *pdbload0(const char *fname, int verbose);
pdbaabb_t *pdbgetaabb(pdbmodel_t *m, int verbose);
#define pdbmdlfree(m) { free(m->at); free(m); }
#define pdbaabbfree(b) { free(b->res); free(b); }

int pdbaaidx(const char *res);
const char *pdbaaname(int i);

#endif


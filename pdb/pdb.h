#include "rv3.h"
#include "util.h"
#ifndef PDB_H__
#define PDB_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct {
  int aid; /* atom index */
  int rid; /* residue index */
  int insert;
  char atnm[8];
  char resnm[8];
  real *x; /* pointer to pdbmodel_t.x */
} pdbatom_t; /* a single atom entry */

typedef struct {
  int natm; /* # of lines == # of atoms */
  int nalloc;
  int nres;
  pdbatom_t *atm; /* array of n or nalloc */
  real (*x)[3];  /* coordinate array of n or nalloc */
  const char *file;
} pdbmodel_t; /* raw data in a pdb model */

/* indices for atoms in standard amino acids */
#define AA_CA     0
#define AA_N      1
#define AA_H2     2
#define AA_H3     3
#define AA_C      4
#define AA_O      5
#define AA_OC1    5
#define AA_OC2    6
#define AA_OXT    6
#define AA_H      7
#define AA_H1     7
#define AA_CB     8
#define AA_HA1    8
#define AA_HA     9
#define AA_HA2    9
#define AA_HA3    10
#define AA_HB     10
#define AA_HB1    10
#define AA_CG2    11
#define AA_HB2    11
#define AA_CG     12
#define AA_HB3    12
#define AA_HG21   12
#define AA_OG     12
#define AA_SG     12
#define AA_HD1    13
#define AA_HG     13
#define AA_HG13   13
#define AA_ND1    13
#define AA_OD1    13
#define AA_OE1    13
#define AA_OG1    13
#define AA_SD     13
#define AA_HD2    14
#define AA_HD21   14
#define AA_HE21   14
#define AA_HE3    14
#define AA_OD2    14
#define AA_OE2    14
#define AA_CD1    15
#define AA_CG1    15
#define AA_HG1    15
#define AA_ND2    15
#define AA_CD     16
#define AA_CE1    16
#define AA_HD22   16
#define AA_NE1    16
#define AA_CE2    17
#define AA_HD11   17
#define AA_HG2    17
#define AA_HG22   17
#define AA_HD12   18
#define AA_HE1    18
#define AA_HE22   18
#define AA_HG23   18
#define AA_NE     18
#define AA_CD2    19
#define AA_HG11   19
#define AA_HG3    19
#define AA_CZ2    20
#define AA_HD13   20
#define AA_HE     20
#define AA_HE2    20
#define AA_HG12   20
#define AA_CE     21
#define AA_CH2    21
#define AA_CZ     21
#define AA_HD23   21
#define AA_HD3    21
#define AA_NE2    21
#define AA_HZ     22
#define AA_HZ2    22
#define AA_NH1    22
#define AA_OH     22
#define AA_HH     23
#define AA_HH11   23
#define AA_HZ3    23
#define AA_HH12   24
#define AA_HH2    24
#define AA_NZ     24
#define AA_CZ3    25
#define AA_HZ1    25
#define AA_NH2    25
#define AA_CE3    26
#define AA_HH21   26
#define AA_HH22   27

#define AA_MAXBIT 28

typedef struct {
  int aa;  /* index of amino acid [0, 20) */
  int nat; /* number of atoms */
  int id[AA_MAXBIT]; /* indices to the coordinate array */
  unsigned flags;
  real *xca, *xn, *xc;
} pdbaar_t; /* amino-acid residues */

typedef struct {
  int nres;
  int natm;
  pdbaar_t *res; /* array of nres */
  real (*x)[3]; /* array of natom */
  const char *file; /* input file */
} pdbaac_t; /* amino-acid chain */

/* generic pdb model */
pdbmodel_t *pdbm_read(const char *fname, int verbose);
int pdbm_write(pdbmodel_t *m, const char *fn);
#define pdbm_free(m) { free(m->atm); free(m->x); free(m); }

enum { PDB_CONTACT_CA, PDB_CONTACT_HEAVY, PDB_CONTACT_ALL }; /* ways of searching contacts */
int *pdbm_contact(pdbmodel_t *pm, double rc, int level, int nearby, int dbg);

/* protein pdb */
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose);
#define pdbaac_free(c) { free(c->res); free(c->x); free(c); }
#define pdbaac_x(c, i, nm) c->x[ c->res[i].id[AA_ ## nm] ]

static const char pdb_aanames_[20][4] = {
  "GLY", "ALA", "VAL", "LEU", "ILE", "PRO",
  "SER", "THR", "CYS", "MET",
  "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", 
  "HIS", "PHE", "TYR", "TRP"};

INLINE int pdbaaidx(const char *res)
{
  int i;
  for (i = 0; i < 20; i++)
    if (strcmp(res, pdb_aanames_[i]) == 0) 
      return i;
  return -1;
}

INLINE const char *pdbaaname(int i)
{
  die_if (i < 0 || i >= 20, "invalid amino acid id %d\n", i);
  return pdb_aanames_[i];
}

/* format atom name, out could be equal to atnm 
 * style 0: atom name first, e.g., "HE21", "CB", or style 1: "1HE2" or " CB " */
INLINE char *pdbm_fmtatom(char *out, char *atnm, int style)
{
  size_t n = strlen(atnm), i;
  char c;

  if (n > 4) { /* a mistake */
    fprintf(stderr, "atom name [%s] too long!\n", atnm);
    exit(1);
  }
  if (style == 0) {
    if (isdigit(c = atnm[0])) { /* rotate the string, such that 1HE2 --> HE21; */
      for (i = 1; i < n; i++) out[i-1] = atnm[i];
      out[n-1] = c;
      out[n] = '\0';
    } else if (out != atnm) strcpy(out, atnm);
  } else if (style == 1) {
    if (n == 4) {
      if (isalpha(atnm[0]) && isdigit(c = atnm[n-1])) { /* HE21 --> 1HE2 */
        for (i = n-1; i > 0; i--) out[i] = atnm[i-1];
        out[0] = c;
        out[n] = '\0';
      } else if (out != atnm) strcpy(out, atnm);
    } else { /* n < 4 */
      if (atnm[0] != ' ') {
        for (i = n; i > 0; i--) out[i] = atnm[i-1];
        out[0] = ' ';
      } else if (out != atnm) strcpy(out, atnm);
      for (i = n+1; i < 4; i++) out[i] = ' ';
      out[4] = '\0';
    }
  }
  return out;
}

#endif


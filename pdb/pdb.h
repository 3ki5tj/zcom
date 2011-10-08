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
enum { PDB_CONTACT_CA, PDB_CONTACT_HEAVY, PDB_CONTACT_ALL };
int *pdbm_contact(pdbmodel_t *pm, double rc, int level, int nearby, int dbg);

/* protein pdb */
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose);
#define pdbaac_free(c) { free(c->res); free(c); }

static const char pdb_aanames_[21][4] = {"---",
  "GLY", "ALA", "VAL", "LEU", "ILE", "PRO",
  "SER", "THR", "CYS", "MET", "PHE", "TYR", "TRP", 
  "ASN", "GLN", "ASP", "GLU", "LYS", "HIS", "ARG"};

INLINE int pdbaaidx(const char *res)
{
  int i;
  for (i = 1; i <= 20; i++)
    if (strcmp(res, pdb_aanames_[i]) == 0) 
      return i;
  return -1;
}

INLINE const char *pdbaaname(int i)
{
  if (i <= 0 || i > 20) {
    fprintf(stderr, "invalid pdb id %d\n", i);
    exit(1);
  }
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


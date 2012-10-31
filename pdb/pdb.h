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
  int insert; /* PDB insertion code, column 27 */
  char atnm[8];
  char resnm[8];
  char elem[4];
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


typedef struct {
  int aa;  /* index of amino acid [0, 20) */
  int nat; /* number of atoms */
  int id[32]; /* indices to the coordinate array */
  unsigned long flags;
  real *xca, *xn, *xc, *xo;
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
INLINE int pdbm_write(pdbmodel_t *m, const char *fn);
#define pdbm_free(m) { free(m->atm); free(m->x); free(m); }

enum { PDB_CONTACT_CA, PDB_CONTACT_HEAVY, PDB_CONTACT_ALL }; /* ways of searching contacts */
int *pdbm_contact(pdbmodel_t *pm, double rc, int level, int nearby, int dbg);

/* protein pdb */
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose);
#define pdbaac_free(c) { free(c->res); free(c->x); free(c); }
#define pdbaac_x(c, i, nm) pdbaac_getx(c, i, #nm)

#define AA_CA         0x1
#define AA_N          0x2
#define AA_C          0x4
#define AA_BACKBONE   (AA_CA | AA_N | AA_C)
#define AA_MAINCHAIN  (AA_BACKBONE | AA_O)
#define AA_H2   29
#define AA_H3   30
#define AA_OXT  31

/* don't edit data in the structure, written by mkdb.py */
struct tag_pdb_aadb {
  const char *resnm;      /* residue name */
  const char *atom[25];   /* atoms */
  const char *sub[11];    /* substitutions */
  unsigned long hvflags;  /* backbone and heavy atom flags */
} pdb_aadb[20] = {
{"GLY", {"CA", "N", "C", "O", "HA1", "HA2", "H", NULL}, {"HA3", "HA1", NULL}, 0xful},
{"ALA", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "HB3", "HA", "H", NULL}, {NULL}, 0x1ful},
{"VAL", {"CA", "N", "C", "O", "CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23", "HA", "H", NULL}, {NULL}, 0x45ful},
{"LEU", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x229ful},
{"ILE", {"CA", "N", "C", "O", "CB", "HB", "CG2", "HG21", "HG22", "HG23", "CG1", "HG11", "HG12", "CD", "HD1", "HD2", "HD3", "HA", "H", NULL}, {"HG13", "HG11", "CD1", "CD", "HD11", "HD1", "HD12", "HD2", "HD13", "HD3", NULL}, 0x245ful},
{"PRO", {"CA", "N", "C", "O", "CD", "HD1", "HD2", "CG", "HG1", "HG2", "CB", "HB1", "HB2", "HA", NULL}, {"HB3", "HB1", "HG3", "HG1", "HD3", "HD1", NULL}, 0x49ful},
{"SER", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "OG", "HG", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x9ful},
{"THR", {"CA", "N", "C", "O", "CB", "HB", "CG2", "HG21", "HG22", "HG23", "OG1", "HG1", "HA", "H", NULL}, {NULL}, 0x45ful},
{"CYS", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "SG", "HG", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x9ful},
{"MET", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "SD", "CE", "HE1", "HE2", "HE3", "HA", "H", NULL}, {"HB3", "HB1", "HG3", "HG1", NULL}, 0xc9ful},
{"ASN", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "OD1", "ND2", "HD21", "HD22", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x39ful},
{"GLN", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "OE1", "NE2", "HE21", "HE22", "HA", "H", NULL}, {"HB3", "HB1", "HG3", "HG1", NULL}, 0x1c9ful},
{"ASP", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "OD1", "OD2", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x39ful},
{"GLU", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "OE1", "OE2", "HA", "H", NULL}, {"HB3", "HB1", "HG3", "HG1", NULL}, 0x1c9ful},
{"LYS", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3", "HA", "H", NULL}, {"HB3", "HB1", "HG3", "HG1", "HD3", "HD1", "HE3", "HE1", NULL}, 0x1249ful},
{"ARG", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22", "HA", "H", NULL}, {"HB3", "HB1", "HG3", "HG1", "HD3", "HD1", NULL}, 0x9a49ful},
{"HIS", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x559ful},
{"PHE", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "HZ", "CE2", "HE2", "CD2", "HD2", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x1559ful},
{"TYR", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "OH", "HH", "CE2", "HE2", "CD2", "HD2", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x2b59ful},
{"TRP", {"CA", "N", "C", "O", "CB", "HB1", "HB2", "CG", "CD1", "HD1", "NE1", "HE1", "CE2", "CZ2", "HZ2", "CH2", "HH2", "CZ3", "HZ3", "CE3", "HE3", "CD2", "HA", "H", NULL}, {"HB3", "HB1", NULL}, 0x2ab59ful}};


INLINE int pdbaaidx(const char *res)
{
  int i;
  if (strcmp(res, "HID") == 0 || strcmp(res, "HIE") == 0 || strcmp(res, "HIP") == 0)
    res = "HIS"; /* quick fix for HIS */
  for (i = 0; i < 20; i++)
    if (strcmp(res, pdb_aadb[i].resnm) == 0) 
      return i;
  return -1;
}

INLINE const char *pdbaaname(int i)
{
  die_if (i < 0 || i >= 20, "invalid amino acid id %d\n", i);
  return pdb_aadb[i].resnm;
}

/* return the index of an atom from */
INLINE int pdbaagetaid(int aa, const char *atnm)
{
  int k;
  for (k = 0; pdb_aadb[aa].atom[k]; k++)
    if (strcmp(pdb_aadb[aa].atom[k], atnm) == 0) return k;
  return -1;
}

/* return the global atom index */
INLINE int pdbaar_getaid(pdbaar_t *r, const char *atnm)
  { int topid = pdbaagetaid(r->aa, atnm); return (topid < 0) ? -1 : r->id[topid]; }
INLINE int pdbaac_getaid(pdbaac_t *c, int i, const char *atnm)
  { return pdbaar_getaid(c->res + i, atnm); }
/* return coordinates */
INLINE real *pdbaac_getx(pdbaac_t *c, int i, const char *atnm)
  { int id = pdbaac_getaid(c, i, atnm); return (id < 0) ? NULL : c->x[id]; }

/* format atom name, out could be equal to atnm 
 * style 0: atom name first, e.g., "HE21", "CB", or style 1: "1HE2" or " CB " */
INLINE char *pdbm_fmtatom(char *out, const char *inp, int style)
{
  size_t n, i;
  char c, cn, atnm[5];
  const char *p;

  die_if (style > 2 || style < 0, "bad format style %d\n", style);
  
  /* copy inp to a buffer without space */
  for (p = inp; *p == ' '; p++) ;
  for (n = 0; n < 4 && p[n] && p[n] != ' '; n++) atnm[n] = p[n];
  atnm[n] = '\0';
  die_if (n == 4 && p[n] && p[n] != ' ', "bad input atom name [%s]\n", atnm);
  
  if (style <= 1) { /* style 0: "H", "CA", "HE21", style 1: " H", " CA", "HE21" */
    if (n == 4) {
      c = atnm[0];
      if (isdigit(c)) { /* rotate the string, such that 1HE2 --> HE21; */
        for (i = 1; i < n; i++) out[i-1] = atnm[i];
        out[n-1] = c;
        out[n] = '\0';
      } else strcpy(out, atnm);
    } else { /* n <= 3 */
      if (style == 0) strcpy(out, atnm);
      else { out[0] = ' '; strcpy(out+1, atnm); }
    }
  } else if (style == 2) { /* style 2: " H", " CA", "1HE2" */
    if (n == 4) {
      c = atnm[0];
      cn = atnm[n - 1];
      if (isalpha(c) && isdigit(cn)) { /* HE21 --> 1HE2 */
        for (i = n-1; i > 0; i--) out[i] = atnm[i-1];
        out[0] = cn;
        out[n] = '\0';
      } else strcpy(out, atnm);
    } else { /* n <= 3 */
      out[0] = ' '; strcpy(out+1, atnm);
    }
  }
  return out;
}

#endif


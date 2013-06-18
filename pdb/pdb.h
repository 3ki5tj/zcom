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
  int unitnm; /* unit is nanometer instead of angstrom */
  const char *file;  /* name of the source file */
} pdbmodel_t; /* raw data in a pdb model */

/* major functions for generic pdb model */
#define pdbm_open(fn, v) pdbm_read(fn, v)
pdbmodel_t *pdbm_read(const char *fn, int verbose);
INLINE int pdbm_write(pdbmodel_t *m, const char *fn);



typedef struct {
  int iaa;  /* index of amino acid [0, 20) */
  int nat; /* number of atoms */
  int id[32]; /* indices to the coordinate array */
  unsigned long flags;
  real *xca, *xn, *xc, *xo; /* pointers to pdbacc->x */
} pdbaar_t; /* amino-acid residues */

typedef struct {
  int nres;
  int natm;
  pdbaar_t *res; /* array of nres */
  real (*x)[3]; /* array of natom */
  int unitnm; /* unit is nanometer instead of angstrom */
  const char *file; /* input file */
} pdbaac_t; /* amino-acid chain */



/* protein pdb */
INLINE pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose);
INLINE pdbaac_t *pdbaac_open(const char *fn, int verbose);
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
{ int topid = pdbaagetaid(r->iaa, atnm); return (topid < 0) ? -1 : r->id[topid]; }

INLINE int pdbaac_getaid(pdbaac_t *c, int i, const char *atnm)
{ return pdbaar_getaid(c->res + i, atnm); }

/* return coordinates */
INLINE real *pdbaac_getx(pdbaac_t *c, int i, const char *atnm)
{ int id = pdbaac_getaid(c, i, atnm); return (id < 0) ? NULL : c->x[id]; }



/* format atom name, `out' can be same as `inp',
 * `style' can be
 * 0: atom name first:   "H",  "CB", "HE21"
 * 1: atom name first:  " H", " CB", "HE21"
 * 2: atom index first: " H", " CB", "1HE2" */
INLINE char *pdbm_fmtatom(char *out, const char *inp, int style)
{
  size_t n, i;
  char c, cn, atnm[5];
  const char *p;

  die_if (style > 2 || style < 0, "bad format style %d\n", style);

  /* copy `inp' to a buffer without space */
  for (p = inp; *p == ' '; p++) ; /* skip the leading spaces */
  for (n = 0; n < 4 && p[n] && p[n] != ' '; n++)
    atnm[n] = p[n]; /* copy without trailing spaces */
  atnm[n] = '\0';
  die_if (n == 4 && p[n] && p[n] != ' ',
      "bad input atom name [%s]\n", atnm);

  if (style <= 1) { /* style 0:  "H",  "CA", "HE21"
                       style 1: " H", " CA", "HE21" */
    if (n == 4) { /* four characters */
      c = atnm[0];
      if (isdigit(c)) { /* rotate the string: 1HE2 --> HE21; */
        for (i = 1; i < n; i++) out[i-1] = atnm[i];
        out[n-1] = c;
        out[n] = '\0';
      } else strcpy(out, atnm);
    } else { /* one to three characters */
      if (style == 0) {
        strcpy(out, atnm);
      } else { /* leave the first character blank */
        out[0] = ' ';
        strcpy(out+1, atnm);
      }
    }
  } else if (style == 2) { /* style 2: " H", " CA", "1HE2" */
    if (n == 4) { /* four characters */
      c = atnm[0];
      cn = atnm[n - 1];
      if (isalpha(c) && isdigit(cn)) { /* HE21 --> 1HE2 */
        for (i = n-1; i > 0; i--) out[i] = atnm[i-1];
        out[0] = cn;
        out[n] = '\0';
      } else strcpy(out, atnm);
    } else { /* n <= 3, leave the first character blank */
      out[0] = ' ';
      strcpy(out+1, atnm);
    }
  }
  return out;
}



/* switch the unit between angstrom and nanometer */
#define pdbm_a2nm(m)    pdbxswitchunit(m->x, m->natm, &m->unitnm, 0)
#define pdbm_nm2a(m)    pdbxswitchunit(m->x, m->natm, &m->unitnm, 1)
#define pdbaac_a2nm(c)  pdbxswitchunit(c->x, c->natm, &c->unitnm, 0)
#define pdbaac_nm2a(c)  pdbxswitchunit(c->x, c->natm, &c->unitnm, 1)

INLINE void pdbxswitchunit(rv3_t *x, int n, int *unitnm, int nm2a)
{
  die_if (nm2a != *unitnm,
    "unit corruption, nm2a %d, unitnm %d\n", nm2a, *unitnm);
  if (x != NULL) {
    real fac = (real)(nm2a ? 10.0 : 0.1);
    int i;

    for (i = 0; i < n; i++) rv3_smul(x[i], fac);
  }
  *unitnm = !(*unitnm); /* switch the unit */
}



/* read raw atom data from pdb */
static pdbmodel_t *pdbm_readpdb(const char *fn)
{
  const int BSIZ = 256;
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *atm;
  int i, ir;
  char s[256], resnm[8] = "";
  float x[3];

  xfopen(fp, fn, "r", return NULL);
  xnew(m, 1);
  m->natm = 0;
  m->nalloc = BSIZ;
  m->file = fn;
  xnew(m->atm, m->nalloc);
  xnew(m->x,   m->nalloc);
  /* read through pdb */
  ir = -1;
  while (fgets(s, sizeof s, fp)) {
    if (strncmp(s, "TER", 3) == 0 ||
        strncmp(s, "ENDMDL", 6) == 0 ||
        strncmp(s, "END", 3) == 0)
      break;
    if (strncmp(s, "ATOM ", 5) != 0)
      continue;
    if (s[16] != ' ' && s[16] != 'A') /* discard alternative position */
      continue;
    i = m->natm;
    if (i >= m->nalloc) {
      m->nalloc += BSIZ;
      xrenew(m->atm, m->nalloc);
      xrenew(m->x,   m->nalloc);
    }
    atm = m->atm + i;
    atm->aid = i;
    atm->insert = s[26];

    /* atom name */
    strip( substr(atm->atnm, s, 12, 4) );
    pdbm_fmtatom(atm->atnm, atm->atnm, 0);
    /* residue name */
    strip( substr(atm->resnm, s, 17, 3) );
    /* residue number */
    sscanf(s+22, "%d", &(atm->rid));
    if (ir == atm->rid && resnm[0] && strcmp(atm->resnm, resnm) != 0) {
      fprintf(stderr, "atom %d, %s, residue %d conflicts %s --> %s, file %s\n",
          i, atm->atnm, ir, resnm, atm->resnm, fn);
    }
    strcpy(resnm, atm->resnm);
    ir = atm->rid;
    /* coordinates */
    sscanf(s+30, "%f%f%f", x, x+1, x+2);
    rv3_make(m->x[i], x[0], x[1], x[2]);
    if (m->unitnm) rv3_smul(m->x[i], (real) 0.1);
    /* element name */
    atm->elem[0] = '\0';
    if (strlen(s) >= 78)
      strip( substr(atm->elem, s, 76, 2) );
    if (atm->elem[0] == '\0') { /* guess */
      atm->elem[0] = atm->atnm[0];
      atm->elem[1] = '\0';
    }
    m->natm++;
  }
  for (i = 0; i < m->natm; i++) /* set atom x */
    m->atm[i].x = m->x[i];
  fclose(fp);
  return m;
}



/* read from GROMACS .gro format */
static pdbmodel_t *pdbm_readgro(const char *fn)
{
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *atm;
  int i, ir;
  char s[256], resnm[8] = "";
  float x[3];

  xfopen(fp, fn, "r", return NULL);
  xnew(m, 1);
  if (fgets(s, sizeof s, fp) == NULL) { /* title line */
    fprintf(stderr, "cannot read the first line of %s\n", fn);
    goto ERR;
  }
  if (fgets(s, sizeof s, fp) == NULL) { /* number of particles */
    fprintf(stderr, "cannot read the second line of %s\n", fn);
    goto ERR;
  }
  if (1 != sscanf(s, "%d", &m->natm)) {
    fprintf(stderr, "no # of atoms, %s\n", fn);
    goto ERR;
  }
  m->nalloc = m->natm;
  m->file = fn;
  xnew(m->atm, m->nalloc);
  xnew(m->x,   m->nalloc);
  ir = -1;
  for (i = 0; i < m->natm; i++) {
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "unable to read atom %d\n", i+1);
      goto ERR;
    }
    atm = m->atm + i;
    atm->aid = i;

    /* atom name */
    strip( substr(atm->atnm, s, 10, 5) );
    pdbm_fmtatom(atm->atnm, atm->atnm, 0);
    /* residue name */
    strip( substr(atm->resnm, s, 5, 5) );
    /* residue number */
    sscanf(s, "%d", &atm->rid);
    if (ir == atm->rid && resnm[0] && strcmp(atm->resnm, resnm) != 0) {
      fprintf(stderr, "atom %d, %s, residue %d conflicts %s --> %s, file %s\n",
          i, atm->atnm, ir, resnm, atm->resnm, fn);
    }
    strcpy(resnm, atm->resnm);
    ir = atm->rid;
    /* coordinates */
    sscanf(s+20, "%f%f%f", x, x+1, x+2);
    rv3_make(m->x[i], x[0], x[1], x[2]);
    if (!m->unitnm) /* x 10 if the unit is nm */
      rv3_smul(m->x[i], (real) 10.);
    atm->x = m->x[i];
    /* guess element name */
    atm->elem[0] = atm->atnm[0];
    atm->elem[1] = '\0';
  }
  fclose(fp);
  return m;
ERR:
  free(m);
  fclose(fp);
  return NULL;
}



/* read pdb */
pdbmodel_t *pdbm_read(const char *fn, int verbose)
{
  int i, j, ir, iro;
  pdbmodel_t *m;
  const char *p;

  p = strrchr(fn, '.');
  if (p != NULL && strcmp(p + 1, "gro") == 0) {
    m = pdbm_readgro(fn);
  } else {
    m = pdbm_readpdb(fn);
  }
  if (m == NULL) return NULL;

  if (verbose)
    printf("%s has %d residues\n", fn, m->atm[m->natm-1].rid);

  /* sort residue indices */
  for (ir = 0, i = 0; i < m->natm; ir++) {
    iro = m->atm[i].rid;
    for (j = i; j < m->natm && m->atm[j].rid == iro &&
        strcmp(m->atm[j].resnm, m->atm[i].resnm) == 0; j++) {
      m->atm[j].rid = ir;
    }
    if (verbose >= 2)
      printf("atoms %d to %d --> residue %s, %d (%d)\n",
        i+1, j, m->atm[i].resnm, iro, ir+1);
    i = j;
  }
  m->nres = ir;

  if (verbose >= 3) /* dump the PDB */
    for (i = 0; i < m->natm; i++) {
      pdbatom_t *atm = m->atm + i;
      printf("%4d %4s %4d %4s %8.3f %8.3f %8.3f\n",
          atm->aid+1, atm->atnm, atm->rid+1, atm->resnm,
          atm->x[0], atm->x[1], atm->x[2]);
    }
  return m;
}


#define pdbm_free(m) pdbm_close(m)

void pdbm_close(pdbmodel_t *m)
{
  if (m->atm) free(m->atm);
  if (m->x) free(m->x);
  free(m);
}


/* write data to file fn */
INLINE int pdbm_write(pdbmodel_t *m, const char *fn)
{
  int i, aid, ATMFMT = 1;
  char atnm[8] = "";
  FILE *fp;
  pdbatom_t *atm;
  real x[3];

  if (m->natm <= 0) return 1;
  xfopen(fp, fn, "w", return 1);
  for (aid = 1, i = 0; i < m->natm; i++) {
    atm = m->atm + i;
    pdbm_fmtatom(atnm, atm->atnm, ATMFMT);
    rv3_copy(x, m->x[i]);
    if (m->unitnm) rv3_smul(x, (real) 10.);
    fprintf(fp, "ATOM  %5d %-4s %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
        aid++, atnm, atm->resnm, atm->rid+1, x[0], x[1], x[2], atm->elem);
  }
  fprintf(fp, "TER   %5d      %-4sA%4d%54s\n", m->natm+1, atm->resnm, atm->rid+1, "");
  fprintf(fp, "END%77s\n", " ");
  fclose(fp);
  return 0;
}



enum { PDB_CONTACT_CA, PDB_CONTACT_HEAVY, PDB_CONTACT_ALL }; /* ways of searching contacts */

INLINE int iscontactatom(int level, const char *atnm)
{
  if (level == PDB_CONTACT_ALL) return 1;
  else if (level == PDB_CONTACT_HEAVY) return atnm[0] != 'H';
  else return strcmp(atnm, "CA") == 0; /* PDB_CONTACT_CA */
}

/* return a nres x nres matrix that defines if two residues are contacts
 * a pair is considered as a contact only if the distance of two
 * specific-type atoms from the two residues is less than rc
 * 'level' : types of atoms to be used in defining atoms
 *           PDB_CONTACT_CA:      only alpha carbon atoms
 *           PDB_CONTACT_HEAVY:   non-hydrogen atoms
 *           PDB_CONTACT_ALL:     include hydrogen atoms
 * 'nearby': # of adjacent resdiues to be excluded from the list
 * */
int *pdbm_contact(pdbmodel_t *pm, double rc, int level, int nearby, int dbg)
{
  int ir, jr, i, j, im, jm, ica, jca, n = pm->natm, nres = pm->nres, ct, cnt = 0;
  pdbatom_t *atm = pm->atm;
  real d, dmin, dca;
  int *iscont;

  xnew(iscont, nres*nres);
  for (ir = 0; ir < nres; ir++) {
    for (jr = ir + nearby; jr < nres; jr++) {
      /* compute the minimal distance between atoms
       * in residues `ir' and `jr' of certain types */
      dmin = 1e9; dca = 0; im = jm = -1;
      /* loop through atoms to collect those of residue id `ir'
       * the loop is inefficient, but this function is rarely called */
      for (i = 0; i < n; i++) {
        if (atm[i].rid != ir) continue;
        if (!iscontactatom(level, atm[i].atnm)) continue;
        ica = iscontactatom(PDB_CONTACT_CA, atm[i].atnm);
        /* loop through atoms to collect those of residue id `jr' */
        for (j = 0; j < n; j++) {
          if (atm[j].rid != jr) continue;
          if (!iscontactatom(level, atm[j].atnm)) continue;
          jca = iscontactatom(PDB_CONTACT_CA, atm[j].atnm);
          d = rv3_dist(atm[i].x, atm[j].x);
          if (d < dmin) { dmin = d; im = i; jm = j; }
          if (ica && jca) dca = d; /* CA distance */
        }
      }
      iscont[ir*nres+jr] = iscont[jr*nres+ir] = ct = (dmin < rc) ? 1 : 0;
      if (ct) cnt++;
      if (dbg && ct) /* print decision */
        printf("[%3d] %s%-3d and %s%-3d: dca %6.3fA dmin %6.3fA (%s:%d, %s:%d)\n",
          cnt, atm[im].resnm, ir+1, atm[jm].resnm, jr+1, dca, dmin,
          atm[im].atnm, im+1, atm[jm].atnm, jm+1);
    }
  }
  return iscont;
}



/* build a `pdbacc_t' structure of amino-acid chain information
 * by parsing `m', which is `pdbmodel_t'
 * the coordinates are duplicated, and `m' can be freed afterwards */
INLINE pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose)
{
  pdbatom_t *atm;
  pdbaac_t *c;
  pdbaar_t *r;
  int i, k, match;
  unsigned long hvflags;

  xnew(c, 1);
  c->nres = m->nres;
  c->natm = m->natm;
  xnew(c->res, c->nres);
  xnew(c->x, m->natm);
  memcpy(c->x, m->x, sizeof(*(c->x)) * m->natm); /* copy coordinates */
  c->unitnm = m->unitnm;
  c->file = m->file;

  for (i = 0; i < m->natm; i++) {
    atm = m->atm + i;
    r = c->res + atm->rid;
    r->iaa = pdbaaidx(atm->resnm);
    if (r->iaa < 0) {
      fprintf(stderr, "unknown amino acid residue %d/%d[%s]\n",
          atm->rid, m->nres, atm->resnm);
      goto ERR;
    }

    /* match index */
    match = 0;
    for (k = 0; pdb_aadb[r->iaa].atom[k] != NULL; k++) {
      if (strcmp(atm->atnm, pdb_aadb[r->iaa].atom[k]) == 0) {
        r->id[k] = i;
        r->flags |= 1ul << k;
        match = 1;
        break;
      }
    }
    if (!match) { /* check terminals */

#define AAMAPIT_(myif, str, nm) \
    myif (strcmp(atm->atnm, str) == 0) { \
      int aid = pdbaagetaid(r->iaa, #nm); \
      r->id[aid] = i; \
      r->flags |= 1ul << (unsigned long) aid; \
      match = 1; }

#define AAMATCH_(myif, str, nm) \
    myif (strcmp(atm->atnm, str) == 0) { \
      r->id[AA_ ## nm] = i; \
      r->flags |= 1ul << AA_##nm; \
      match = 1; }

      AAMAPIT_(if, "H1", H)
      AAMATCH_(else if, "H2",   H2)
      AAMATCH_(else if, "H3",   H3)
      AAMATCH_(else if, "OXT",  OXT)
      AAMAPIT_(else if, "OC1",  O)
      AAMATCH_(else if, "OC2",  OXT)
      AAMAPIT_(else if, "O1",   O)
      AAMATCH_(else if, "O2",   OXT)
      else { /* check substitutions */
        for (k = 0; pdb_aadb[r->iaa].sub[k] != NULL; k += 2)
          if (strcmp(atm->atnm, pdb_aadb[r->iaa].sub[k]) == 0) {
            int aid = pdbaagetaid(r->iaa, pdb_aadb[r->iaa].sub[k+1]);
            r->id[aid] = i;
            r->flags |= 1ul << (unsigned) aid;
            match = 1;
            break;
          }
      }
    }
    if (!match)
      printf("unknown atom %s:%d res %s%d\n",
          atm->atnm, i+1, atm->resnm, atm->rid+1);
  }

#define pdbaac_pmiss_(xflags) { \
  unsigned long miss = (r->flags ^ xflags) & xflags; \
  fprintf(stderr, "file %s, residue %s%d misses atom(s): ", \
          c->file, pdbaaname(r->iaa), i+1); \
  for (k = 0; pdb_aadb[r->iaa].atom[k] != NULL; k++) { \
    if (miss & (1ul << k)) \
      fprintf(stderr, "%s ", pdb_aadb[r->iaa].atom[k]); } \
  fprintf(stderr, "\n"); }

  /* checking integrity */
  for (i = 0; i < c->nres; i++) {
    r = c->res + i;
    hvflags = pdb_aadb[r->iaa].hvflags;
    if ((r->flags & AA_BACKBONE) != AA_BACKBONE) {
      pdbaac_pmiss_(AA_BACKBONE);
      goto ERR;
    } else if ((r->flags & hvflags) != hvflags) {
      pdbaac_pmiss_(hvflags);
    }
    r->xn  = pdbaac_x(c, i, N);
    r->xca = pdbaac_x(c, i, CA);
    r->xc  = pdbaac_x(c, i, C);
    r->xo  = pdbaac_x(c, i, O);
  }

  /* check bond-length, assume backbone are present */
  for (i = 0; i < c->nres; i++) {
    double x;
    r = c->res + i;
    if (i > 0) {
      x = rv3_dist(pdbaac_x(c, i-1, C), r->xn);
      if (c->unitnm) x *= 10.;
      if (x < .3 || x > 2.3) {
        if (verbose) {
          const char *aap = pdbaaname(c->res[i - 1].iaa);
          const char *aa = pdbaaname(r->iaa);
          fprintf(stderr, "%s: C-N bond between %d (%s) and %d (%s) is broken %g, insert break\n",
            c->file, i, aap, i+1, aa, x);
          goto ERR;
        }
      }
    }
    x = rv3_dist(r->xn, r->xca);
    if (c->unitnm) x *= 10.;
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: N-CA bond of residue %d (%s) is broken %g\n",
          c->file, i+1, pdbaaname(r->iaa), x);
        fprintf(stderr, "N : %8.3f %8.3f %8.3f\n", r->xn[0], r->xn[1], r->xn[2]);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
      }
      goto ERR;
    }
    x = rv3_dist(r->xca, r->xc);
    if (c->unitnm) x *= 10.;
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: CA-C bond of residue %d (%s) is broken %g\n",
          c->file, i+1, pdbaaname(r->iaa), x);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
        fprintf(stderr, "C : %8.3f %8.3f %8.3f\n", r->xc[0], r->xc[1], r->xc[2]);
      }
      goto ERR;
    }
  }
  if (verbose >= 3) {
    for (i = 0; i < c->nres; i++)
      printf("%4d: %s\n", i+1, pdbaaname(c->res[i].iaa));
  }

  return c;
ERR:
  free(c->res);
  free(c->x);
  free(c);
  return NULL;
}



/* create `pdbaac_t' from a PDB file */
INLINE pdbaac_t *pdbaac_open(const char *fn, int verbose)
{
  pdbmodel_t *m;
  pdbaac_t *c;

  m = pdbm_read(fn, verbose);
  c = pdbaac_parse(m, verbose);
  pdbm_free(m);
  return c;
}



#define pdbaac_free(c) pdbaac_close(c)

INLINE void pdbaac_close(pdbaac_t *c)
{
  if (c->res) free(c->res);
  if (c->x) free(c->x);
  free(c);
}


/* parse helices, return the number of helices `nse'
 * (*pse)[0..nse*2 - 1] are start and finishing indices of helices */
INLINE int pdbaac_parsehelices(pdbaac_t *c, int **pse)
{
  int i, nse = 0, is, it, nres = c->nres;
  int aa, aagly, aapro;
  int *se, *ishx, q[5];
  double phi, psi;

  /* A. make an array of nres, identify if each residue is helix */
  xnew(ishx, nres);
  ishx[0] = ishx[nres-1] = 0;
  for (i = 1; i < nres-1; i++) {
    /* make a local 5-tuple */
    die_if ((q[0] = pdbaac_getaid(c, i-1, "C"))  < 0, "no C  of %d\n", i-1);
    die_if ((q[1] = pdbaac_getaid(c, i,   "N"))  < 0, "no N  of %d\n", i);
    die_if ((q[2] = pdbaac_getaid(c, i,   "CA")) < 0, "no CA of %d\n", i);
    die_if ((q[3] = pdbaac_getaid(c, i,   "C"))  < 0, "no C  of %d\n", i);
    die_if ((q[4] = pdbaac_getaid(c, i+1, "N"))  < 0, "no N  of %d\n", i+1);
    phi = rv3_dih(c->x[ q[0] ], c->x[ q[1] ], c->x[ q[2] ], c->x[ q[3] ],
                  NULL, NULL, NULL, NULL);
    psi = rv3_dih(c->x[ q[1] ], c->x[ q[2] ], c->x[ q[3] ], c->x[ q[4] ],
                  NULL, NULL, NULL, NULL);
    ishx[i] = (phi < 0 && psi > -100*M_PI/180 && psi < 80*M_PI/180);
  }

  /* B. searching for segments
   * make 2*pro->ngrp for start/end of each segment
   * range of segment k is se[2*k] <= id < se[2*k+1] */
  xnew(se, 2);
  aagly = pdbaaidx("GLY");
  aapro = pdbaaidx("PRO");
  for (i = 0, is = 0; i < nres; ) { /* try to find the helices */
    while (i < nres && !ishx[i]) i++;
    if (i >= nres) break; /* no more helices */
    is = i;
    while (ishx[i] && i < nres) i++;
    it = i;
    for (; is < it; is++) { /* skip terminal GLY and PRO */
      aa = c->res[is].iaa;
      if (aa != aagly && aa != aapro)  break;
    }
    for (; it > is; it--) {
      aa = c->res[it - 1].iaa;
      if (aa != aagly && aa != aapro) break;
    }
    if (it - is >= 4) { /* successfully find a helical segment */
      xrenew(se, 2*(nse+1));
      se[2*nse] = is;
      se[2*nse+1] = it;
      nse++;
    } else { } /* just let go, don't increment nse */
  }
  free(ishx);
  *pse = se;
  return nse;
}


#endif


#ifndef PDB_C__
#define PDB_C__
#include "pdb.h"

/* read raw atom data from pdb */
pdbmodel_t *pdbm_read(const char *fname, int verbose)
{
  const int BSIZ = 256;
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *atm;
  int i, j, ir, iro;
  char s[256], resnm[8];
  float x[3];

  xfopen(fp, fname, "r", return NULL);
  xnew(m, 1);
  m->natm = 0;
  m->nalloc = BSIZ;
  m->file = fname;
  xnew(m->atm, m->nalloc);
  xnew(m->x,   m->nalloc);
  /* read through pdb */
  /*nline = 0; */
  ir = -1;
  while (fgets(s, sizeof s, fp)) {
    /*nline++; */
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

#define sncpy_(dest, src, n) { \
  char *p; \
  memcpy(dest, src, n); \
  *(dest+n) = '\0'; \
  for (p = dest; isspace(*p); p++) ; \
  if (*p == '\0') *dest = '\0'; \
  else memmove(dest, p, strlen(p)+1); \
  for (p = dest + strlen(dest) - 1; p >= dest && isspace(*p); p--) *p = '\0'; }

    /* atom name */
    sncpy_(atm->atnm, s+12, 4);
    pdbm_fmtatom(atm->atnm, atm->atnm, 0);
    /* residue name */
    sncpy_(atm->resnm, s+17, 3);
#undef sncpy_
    /* residue number */
    sscanf(s+22, "%d", &(atm->rid));
    if (ir == atm->rid && strcmp(atm->resnm, resnm) != 0) {
      /*
      fprintf(stderr, "residue %d conflicts %s --> %s at line %d, file %s\n",
          ir, resnm, atm->resnm, nline, fname);
      */
    }
    strcpy(resnm, atm->resnm);
    ir = atm->rid;
    /* coordinates */
    sscanf(s+30, "%f%f%f", x, x+1, x+2);
    rv3_make(m->x[i], x[0], x[1], x[2]);
    m->natm++;
  }
  for (i = 0; i < m->natm; i++) /* set atom x */
    m->atm[i].x = m->x[i];
  if (verbose)
    printf("%s has %d residues\n", fname, m->atm[m->natm-1].rid);
  /* offset the residue id */
  for (ir = 0, i = 0; i < m->natm; ir++) {
    iro = m->atm[i].rid;
    for (j = i; j < m->natm && m->atm[j].rid == iro &&
        m->atm[j].insert == m->atm[i].insert &&
        strcmp(m->atm[j].resnm, m->atm[i].resnm) == 0; j++) {
      m->atm[j].rid = ir;
    }
    if (verbose >= 2)
      printf("%d to %d is set for residue %s, %d to %d\n",
        i+1, j+1, m->atm[i].resnm, iro, ir+1);
    i = j;
  }
  m->nres = ir;
  if (verbose >= 3) {
    for (i = 0; i < m->natm; i++) {
      atm = m->atm + i;
      printf("%4d %4s %4d %4s %8.3f %8.3f %8.3f\n",
          atm->aid+1, atm->atnm, atm->rid+1, atm->resnm, 
          atm->x[0], atm->x[1], atm->x[2]);
    }
  }
  fclose(fp);
  return m;
}

/* write data to file fn */
int pdbm_write(pdbmodel_t *m, const char *fn)
{
  int i, aid, atp;
  char atnm[8] = "";
  FILE *fp;
  pdbatom_t *atm;
  real *x;

  xfopen(fp, fn, "w", return 1);
  for (aid = 1, i = 0; i < m->natm; i++) {
    atm = m->atm + i;
    atp = atm->atnm[0];
    pdbm_fmtatom(atnm, atm->atnm, 1);
    x = m->x[i];
    fprintf(fp, "ATOM  %5d %4s %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00           %c  \n", 
        aid++, atnm, atm->resnm, atm->rid+1, x[0], x[1], x[2], atp);
  }
  fprintf(fp, "END%77s\n", " ");
  fclose(fp);
  return 0;
}

INLINE int iscontactatom(int level, const char *atnm)
{
  if (level == PDB_CONTACT_ALL) return 1;
  else if (level == PDB_CONTACT_HEAVY) return atnm[0] != 'H';
  else return strcmp(atnm, "CA") == 0; /* PDB_CONTACT_CA */
}

/* return a nres x nres array to indicate if two residues form a contact
 * it is a contact only if the distance of any two atoms from the two residues 
 * is less than rc,
 * 'level' : one of the PDB_CONTACT_XX values above.
 * 'nearby': # of adjacent resdiues to be excluded from the list
 * */
int *pdbm_contact(pdbmodel_t *pm, double rc, int level, int nearby, int dbg)
{
  int ir, jr, i, j, im, jm, ica, jca, n = pm->natm, nres = pm->nres, ct;
  pdbatom_t *atm = pm->atm;
  real d, dmin, dca;
  int *ds;

  xnew(ds, nres*nres);
  for (ir = 0; ir < nres; ir++) {
    for (jr = ir+1; jr < nres; jr++) {
      /* compute the minimal distance between ir and jr */
      dmin = 1e9; dca = 0; im = jm = -1;
      for (i = 0; i < n; i++) {
        if (atm[i].rid != ir) continue;
        if (!iscontactatom(level, atm[i].atnm)) continue;
        ica = iscontactatom(PDB_CONTACT_CA, atm[i].atnm);
        for (j = 0; j < n; j++) {
          if (atm[j].rid != jr) continue;
          if (!iscontactatom(level, atm[j].atnm)) continue;
          jca = iscontactatom(PDB_CONTACT_CA, atm[j].atnm);
          d = rv3_dist(atm[i].x, atm[j].x);
          if (d < dmin) { dmin = d; im = i; jm = j; }
          if (ica && jca) dca = d; /* CA distance */
        }
      }
      ds[ir*nres+jr] = ds[jr*nres+ir] = ct = (dmin < rc) ? 1 : 0;
      if (dbg > 1 || (dbg && ct))/* print decision */
        printf("[%d] %s%-3d and %s%-3d: dca %6.3fA dmin %6.3fA (%s:%d, %s:%d)\n", 
          ct, atm[im].resnm, ir+1, atm[jm].resnm, jr+1, dca, dmin, 
          atm[im].atnm, im+1, atm[jm].atnm, jm+1);
    }
  }

  /* exclude nearby residues */
  for (ir = 0; ir < nres; ir++)
    for (jr = ir+1; jr <= ir+nearby && jr < nres; jr++)
      ds[ir*nres + jr] = ds[jr*nres + ir] = 0;
  return ds;
}

/* get amino acid chain by parsing pdbmodel_t m */
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose)
{
  pdbatom_t *atm;
  pdbaac_t *c;
  pdbaar_t *r;
  int i, k;
  const char *nm;
  static struct { unsigned bit; const char *nm; } aaatomtypes[] = {
  {AA_CA,   "CA"  }, {AA_N,    "N"   }, {AA_H2,   "H2"  }, {AA_H3,   "H3"  }, {AA_C,    "C"   }, 
  {AA_O,    "O"   }, {AA_OC1,  "OC1" }, {AA_OC2,  "OC2" }, {AA_OXT,  "OXT" }, {AA_H,    "H"   }, 
  {AA_H1,   "H1"  }, {AA_CB,   "CB"  }, {AA_HA1,  "HA1" }, {AA_HA,   "HA"  }, {AA_HA2,  "HA2" }, 
  {AA_HA3,  "HA3" }, {AA_HB,   "HB"  }, {AA_HB1,  "HB1" }, {AA_CG2,  "CG2" }, {AA_HB2,  "HB2" }, 
  {AA_CG,   "CG"  }, {AA_HB3,  "HB3" }, {AA_HG21, "HG21"}, {AA_OG,   "OG"  }, {AA_SG,   "SG"  }, 
  {AA_HD1,  "HD1" }, {AA_HG,   "HG"  }, {AA_HG13, "HG13"}, {AA_ND1,  "ND1" }, {AA_OD1,  "OD1" }, 
  {AA_OE1,  "OE1" }, {AA_OG1,  "OG1" }, {AA_SD,   "SD"  }, {AA_HD2,  "HD2" }, {AA_HD21, "HD21"}, 
  {AA_HE21, "HE21"}, {AA_HE3,  "HE3" }, {AA_OD2,  "OD2" }, {AA_OE2,  "OE2" }, {AA_CD2,  "CD2" }, 
  {AA_CG1,  "CG1" }, {AA_HG1,  "HG1" }, {AA_ND2,  "ND2" }, {AA_CD,   "CD"  }, {AA_CD1,  "CD1" }, 
  {AA_CE1,  "CE1" }, {AA_HD22, "HD22"}, {AA_HG2,  "HG2" }, {AA_HG22, "HG22"}, {AA_NE1,  "NE1" }, 
  {AA_CE2,  "CE2" }, {AA_HD11, "HD11"}, {AA_HG23, "HG23"}, {AA_HG3,  "HG3" }, {AA_HD12, "HD12"}, 
  {AA_HE1,  "HE1" }, {AA_HE22, "HE22"}, {AA_HG11, "HG11"}, {AA_NE,   "NE"  }, {AA_CZ2,  "CZ2" }, 
  {AA_HD13, "HD13"}, {AA_HE,   "HE"  }, {AA_HE2,  "HE2" }, {AA_HG12, "HG12"}, {AA_CE,   "CE"  }, 
  {AA_CH2,  "CH2" }, {AA_CZ,   "CZ"  }, {AA_HD23, "HD23"}, {AA_HD3,  "HD3" }, {AA_NE2,  "NE2" }, 
  {AA_HZ,   "HZ"  }, {AA_HZ2,  "HZ2" }, {AA_NH1,  "NH1" }, {AA_OH,   "OH"  }, {AA_HH,   "HH"  }, 
  {AA_HH11, "HH11"}, {AA_HZ3,  "HZ3" }, {AA_HH12, "HH12"}, {AA_HH2,  "HH2" }, {AA_NZ,   "NZ"  }, 
  {AA_CZ3,  "CZ3" }, {AA_HZ1,  "HZ1" }, {AA_NH2,  "NH2" }, {AA_CE3,  "CE3" }, {AA_HH21, "HH21"}, 
  {AA_HH22, "HH22"}, {-1, NULL}};

  xnew(c, 1);
  c->nres = m->nres;
  c->natm = m->natm;
  xnew(c->res, c->nres);
  xnew(c->x, m->natm);
  memcpy(c->x, m->x, sizeof(*(c->x))*m->natm); /* copy coordinates */
  c->file = m->file;

  for (i = 0; i < m->natm; i++) {
    atm = m->atm + i;
    r = c->res + atm->rid;
    r->aa = pdbaaidx(atm->resnm);
    if (r->aa < 0) {
      fprintf(stderr, "unknown amino acid residue %d/%d[%s]\n", 
          atm->rid, m->nres, atm->resnm);
      goto ERR;
    }

    /* match index */
    for (k = 0; (nm = aaatomtypes[k].nm) != NULL; k++) {
      int bit = aaatomtypes[k].bit;
      if (strcmp(atm->atnm, nm) == 0) {
        r->id[bit] = i;
        r->flags |= 1u << bit;
        break;
      }
    }
    if (nm == NULL) {
      printf("Error: unknown atom %s-%d res %s-%d\n", atm->atnm, i, atm->resnm, atm->rid);
    }
  }

  /* checking integrity */
  for (i = 0; i < c->nres; i++) {
    unsigned bbflags = (1 << AA_CA)|(1 << AA_N)|(1 << AA_C);
    r = c->res + i;
    if ((r->flags & bbflags) != bbflags) {
      if (verbose)
        fprintf(stderr, "%s: coordinates of residue %d/%d (%s) is incomplete (0x%X).\n",
          c->file, i+1, c->nres, pdbaaname(r->aa), r->flags);
    }
    r->xn = pdbaac_x(c, i, N);
    r->xca = pdbaac_x(c, i, CA);
    r->xc = pdbaac_x(c, i, C);
  }
  /* check bond-length */
  for (i = 0; i < c->nres; i++) {
    real x;
    r = c->res + i;
    if (i > 0) {
      x = rv3_dist(pdbaac_x(c, i-1, C), r->xn);
      if (x < .3 || x > 2.3) {
        if (verbose) {
          const char *aap = pdbaaname(c->res[i-1].aa), *aa = pdbaaname(r->aa);
          fprintf(stderr, "%s: C-N bond between %d (%s) and %d (%s) is broken %g, insert break\n",
            c->file, i, aap, i+1, aa, x);
          goto ERR;
        }
      }
    }
    x = rv3_dist(r->xn, r->xca);
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: N-CA bond of residue %d (%s) is broken %g\n",
          c->file, i+1, pdbaaname(r->aa), x);
        fprintf(stderr, "N : %8.3f %8.3f %8.3f\n", r->xn[0], r->xn[1], r->xn[2]);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
      }
      goto ERR;
    }
    x = rv3_dist(r->xca, r->xc);
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: CA-C bond of residue %d (%s) is broken %g\n",
          c->file, i+1, pdbaaname(r->aa), x);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
        fprintf(stderr, "C : %8.3f %8.3f %8.3f\n", r->xc[0], r->xc[1], r->xc[2]);
      }
      goto ERR;
    }
  }
  if (verbose >= 3) {
    for (i = 0; i < c->nres; i++)
      printf("%4d: %s\n", i+1, pdbaaname(c->res[i].aa));
  }
  return c;
ERR:
  free(c->res);
  free(c);
  return NULL;
}

#endif


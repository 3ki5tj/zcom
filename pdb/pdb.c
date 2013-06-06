#ifndef PDB_C__
#define PDB_C__
#include "pdb.h"


/* switch the unit between angstrom and nanometer */
#define pdbm_a2nm(m) pdbm_switchunit(m, 0)
#define pdbm_nm2a(m) pdbm_switchunit(m, 1)
INLINE void pdbm_switchunit(pdbmodel_t *m, int nm2a)
{
  die_if (nm2a != m->unitnm,
    "unit corruption, nm2a %d, unitnm %d\n", nm2a, m->unitnm);
  if (m->x != NULL) {
    real fac = (real)(nm2a ? 10.0 : 0.1);
    int i;

    for (i = 0; i < m->natm; i++) rv3_smul(m->x[i], fac);
  }
  m->unitnm = !m->unitnm; /* switch the unit */
}



/* read raw atom data from pdb */
static pdbmodel_t *pdbm_readpdb(const char *fname)
{
  const int BSIZ = 256;
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *atm;
  int i, ir;
  char s[256], resnm[8] = "";
  float x[3];

  xfopen(fp, fname, "r", return NULL);
  xnew(m, 1);
  m->natm = 0;
  m->nalloc = BSIZ;
  m->file = fname;
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
          i, atm->atnm, ir, resnm, atm->resnm, fname);
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
static pdbmodel_t *pdbm_readgro(const char *fname)
{
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *atm;
  int i, ir;
  char s[256], resnm[8] = "";
  float x[3];

  xfopen(fp, fname, "r", return NULL);
  xnew(m, 1);
  if (fgets(s, sizeof s, fp) == NULL) { /* title line */
    fprintf(stderr, "cannot read the first line of %s\n", fname);
    goto ERR;
  }
  if (fgets(s, sizeof s, fp) == NULL) { /* number of particles */
    fprintf(stderr, "cannot read the second line of %s\n", fname);
    goto ERR;
  }
  if (1 != sscanf(s, "%d", &m->natm)) {
    fprintf(stderr, "no # of atoms, %s\n", fname);
    goto ERR;
  }
  m->nalloc = m->natm;
  m->file = fname;
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
          i, atm->atnm, ir, resnm, atm->resnm, fname);
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
pdbmodel_t *pdbm_read(const char *fname, int verbose)
{
  int i, j, ir, iro;
  pdbmodel_t *m;
  const char *p;

  p = strrchr(fname, '.');
  if (p != NULL && strcmp(p + 1, "gro") == 0) {
    m = pdbm_readgro(fname);
  } else {
    m = pdbm_readpdb(fname);
  }
  if (m == NULL) return NULL;

  if (verbose)
    printf("%s has %d residues\n", fname, m->atm[m->natm-1].rid);

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
    for (jr = ir + nearby + 1; jr < nres; jr++) {
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
pdbaac_t *pdbaac_parse(pdbmodel_t *m, int verbose)
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
    r->aa = pdbaaidx(atm->resnm);
    if (r->aa < 0) {
      fprintf(stderr, "unknown amino acid residue %d/%d[%s]\n",
          atm->rid, m->nres, atm->resnm);
      goto ERR;
    }

    /* match index */
    match = 0;
    for (k = 0; pdb_aadb[r->aa].atom[k] != NULL; k++)
      if (strcmp(atm->atnm, pdb_aadb[r->aa].atom[k]) == 0) {
        r->id[k] = i;
        r->flags |= 1ul << k;
        match = 1;
        break;
      }
    if (!match) { /* check terminals */
#define AAMAPIT_(lead, str, nm) lead (strcmp(atm->atnm, str) == 0) \
    { int aid = pdbaagetaid(r->aa, #nm); r->id[aid] = i; r->flags |= 1ul << (unsigned long) aid; match = 1; }
#define AAMATCH_(lead, str, nm) lead (strcmp(atm->atnm, str) == 0) \
    { r->id[AA_ ## nm] = i; r->flags |= 1ul << AA_##nm; match = 1; }
      AAMAPIT_(if, "H1", H)
      AAMATCH_(else if, "H2",   H2)
      AAMATCH_(else if, "H3",   H3)
      AAMATCH_(else if, "OXT",  OXT)
      AAMAPIT_(else if, "OC1",  O)
      AAMATCH_(else if, "OC2",  OXT)
      AAMAPIT_(else if, "O1",   O)
      AAMATCH_(else if, "O2",   OXT)
      else { /* check substitutions */
        for (k = 0; pdb_aadb[r->aa].sub[k] != NULL; k += 2)
          if (strcmp(atm->atnm, pdb_aadb[r->aa].sub[k]) == 0) {
            int aid = pdbaagetaid(r->aa, pdb_aadb[r->aa].sub[k+1]);
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
  fprintf(stderr, "file %s, residue %s%d misses atom(s): ", c->file, pdbaaname(r->aa), i+1); \
  for (k = 0; pdb_aadb[r->aa].atom[k] != NULL; k++) { \
    if (miss & (1ul << k)) fprintf(stderr, "%s ", pdb_aadb[r->aa].atom[k]); } \
  fprintf(stderr, "\n"); }

  /* checking integrity */
  for (i = 0; i < c->nres; i++) {
    r = c->res + i;
    hvflags = pdb_aadb[r->aa].hvflags;
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
          const char *aap = pdbaaname(c->res[i-1].aa), *aa = pdbaaname(r->aa);
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
          c->file, i+1, pdbaaname(r->aa), x);
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
  free(c->x);
  free(c);
  return NULL;
}

/* parse helices, return the number of helices `nse'
 * (*pse)[0..nse*2 - 1] are start and finishing indices of helices */
INLINE int pdbaac_parsehelices(pdbaac_t *c, int **pse)
{
  int i, nse = 0, is, it, nres = c->nres;
  int aa, aagly, aapro;
  int *se, *ishx, quin[5];
  double phi, psi;

  /* A. make an array of nres, identify if each residue is helix */
  xnew(ishx, nres);
  ishx[0] = ishx[nres-1] = 0;
  for (i = 1; i < nres-1; i++) {
    /* make local quintuple */
    quin[0] = pdbaac_getaid(c, i-1, "C");  die_if(quin[0] < 0, "no C  of %d\n", i-1);
    quin[1] = pdbaac_getaid(c, i,   "N");  die_if(quin[1] < 0, "no N  of %d\n", i);
    quin[2] = pdbaac_getaid(c, i,   "CA"); die_if(quin[2] < 0, "no CA of %d\n", i);
    quin[3] = pdbaac_getaid(c, i,   "C");  die_if(quin[3] < 0, "no C  of %d\n", i);
    quin[4] = pdbaac_getaid(c, i+1, "N");  die_if(quin[4] < 0, "no N  of %d\n", i+1);
    phi = rv3_calcdihv(NULL, c->x, quin, 0);
    psi = rv3_calcdihv(NULL, c->x, quin+1, 0);
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
      aa = c->res[is].aa;
      if (aa != aagly && aa != aapro)  break;
    }
    for (; it > is; it--) {
      aa = c->res[it - 1].aa;
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


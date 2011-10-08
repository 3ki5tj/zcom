#include "rv3.h"
#ifndef PDB_C__
#define PDB_C__
#include "pdb.h"

/* read raw atom data from pdb */
pdbmodel_t *pdbm_read(const char *fname, int verbose)
{
  const int bsiz = 256;
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *at;
  int i, j, ir, iro;
  char s[256], resnm[8];
  float x[3];

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fname);
    return NULL;
  }
  if ((m = calloc(1, sizeof(*m))) == NULL) {
    fprintf(stderr, "no memory for pdbmodel\n");
    return NULL;
  }
  m->n = 0;
  m->nalloc = bsiz;
  m->file = fname;
  if ((m->at = calloc(m->nalloc, sizeof(*(m->at)))) == NULL) {
    fprintf(stderr, "no memory for pdbatoms\n");
    free(m);
    return NULL;
  }
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
    i = m->n;
    if (i >= m->nalloc) {
      m->nalloc += bsiz;
      if ((m->at = realloc(m->at, m->nalloc*sizeof(*(m->at)))) == NULL) {
        fprintf(stderr, "no memory to expand m->at");
        free(m);
        return NULL;
      }
    }
    at = m->at + i;
    at->aid = i;
    at->insert = s[26];

#define sncpy_(dest, src, n) { \
  char *p; \
  memcpy(dest, src, n); \
  *(dest+n) = '\0'; \
  for (p = dest; isspace(*p); p++) ; \
  if (*p == '\0') *dest = '\0'; \
  else memmove(dest, p, strlen(p)+1); \
  for (p = dest + strlen(dest) - 1; p >= dest && isspace(*p); p--) *p = '\0'; }

    /* atom name */
    sncpy_(at->atnm, s+12, 4);
    pdbm_fmtatom(at->atnm, at->atnm, 0);
    /* residue name */
    sncpy_(at->resnm, s+17, 3);
#undef sncpy_
    /* residue number */
    sscanf(s+22, "%d", &(at->rid));
    if (ir == at->rid && strcmp(at->resnm, resnm) != 0) {
      /*
      fprintf(stderr, "residue %d conflicts %s --> %s at line %d, file %s\n",
          ir, resnm, at->resnm, nline, fname);
      */
    }
    strcpy(resnm, at->resnm);
    ir = at->rid;
    /* coordinates */
    sscanf(s+30, "%f%f%f", x, x+1, x+2);
    rv3_make(at->x, x[0], x[1], x[2]);
    m->n++;
  }
  if (verbose)
    printf("%s has %d residues\n", fname, m->at[m->n-1].rid);
  /* offset the residue id */
  for (ir = 0, i = 0; i < m->n; ir++) {
    iro = m->at[i].rid;
    for (j = i; j < m->n && m->at[j].rid == iro &&
        m->at[j].insert == m->at[i].insert &&
        strcmp(m->at[j].resnm, m->at[i].resnm) == 0; j++) {
      m->at[j].rid = ir;
    }
    if (verbose >= 2)
      printf("%d to %d is set from residue %s, %d to %d\n",
        i+1, j+1, m->at[i].resnm, iro, ir+1);
    i = j;
  }
  m->nres = ir;
  if (verbose >= 3) {
    for (i = 0; i < m->n; i++) {
      at = m->at + i;
      printf("%4d %4s %4d %4s %8.3f %8.3f %8.3f\n",
          at->aid+1, at->atnm, at->rid+1, at->resnm, 
          at->x[0], at->x[1], at->x[2]);
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
  pdbatom_t *at;

  if ((fp = fopen(fn, "w")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fn);
    return 1;
  }
  for (aid = 1, i = 0; i < m->n; i++) {
    at = m->at + i;
    atp = at->atnm[0];
    pdbm_fmtatom(atnm, at->atnm, 1);
    fprintf(fp, "ATOM  %5d %4s %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00           %c  \n", 
        aid++, atnm, at->resnm, at->rid+1, at->x[0], at->x[1], at->x[2], atp);
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
  int ir, jr, i, j, im, jm, ica, jca, n = pm->n, nres = pm->nres;
  pdbatom_t *at = pm->at;
  real d, dmin, dca;
  int *ds;

  if ((ds = calloc(nres*nres, sizeof(*ds))) == NULL) {
    fprintf(stderr, "cannot allocate contact array\n");
    return NULL;
  }
  for (ir = 0; ir < nres; ir++) {
    for (jr = ir+1; jr < nres; jr++) {
      /* compute the minimal distance between ir and jr */
      dmin = 1e9; dca = 0; im = jm = -1;
      for (i = 0; i < n; i++) {
        if (at[i].rid != ir) continue;
        if (!iscontactatom(level, at[i].atnm)) continue;
        ica = iscontactatom(PDB_CONTACT_CA, at[i].atnm);
        for (j = 0; j < n; j++) {
          if (at[j].rid != jr) continue;
          if (!iscontactatom(level, at[j].atnm)) continue;
          jca = iscontactatom(PDB_CONTACT_CA, at[j].atnm);
          d = rv3_dist(at[i].x, at[j].x);
          if (d < dmin) { dmin = d; im = i; jm = j; }
          if (ica && jca) dca = d; /* CA distance */
        }
      }
      ds[ir*nres+jr] = ds[jr*nres+ir] = (dmin < rc) ? 1 : 0;
      if (dbg && ds[ir*nres + jr]) /* print decision */
        printf("%s%-3d and %s%-3d: dca %6.3fA dmin %6.3fA (%s:%d, %s:%d)\n", 
          at[im].resnm, ir+1, at[jm].resnm, jr+1, dca, dmin, 
          at[im].atnm, im+1, at[jm].atnm, jm+1);
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
  pdbatom_t *at;
  pdbaac_t *c;
  pdbaar_t *r;
  int i;

  if ((c = calloc(1, sizeof(*c))) == NULL) {
    fprintf(stderr, "no memory for aabb\n");
    return NULL;
  }
  c->nres = m->nres;
  if ((c->res = calloc(c->nres, sizeof(*(c->res)))) == NULL) {
    fprintf(stderr, "no memory for residues\n");
    free(c);
    return NULL;
  }
  c->file = m->file;
  for (i = 0; i < m->n; i++) {
    at = m->at + i;
    r = c->res + at->rid;
    r->aa = pdbaaidx(at->resnm);
    if (r->aa < 0) {
      fprintf(stderr, "unknown amino acid residue %d/%d[%s]\n", 
          at->rid, m->nres, at->resnm);
      goto ERR;
    }
    if (strcmp(at->atnm, "N") == 0) {
      memcpy(r->xn, at->x, sizeof(at->x[0])*3);
      r->flags |= HAS_N;
    } else if (strcmp(at->atnm, "CA") == 0) {
      memcpy(r->xca, at->x, sizeof(at->x[0])*3);
      r->flags |= HAS_CA;
    } else if (strcmp(at->atnm, "C") == 0) {
      memcpy(r->xc, at->x, sizeof(at->x[0])*3);
      r->flags |= HAS_C;
    } else if (strcmp(at->atnm, "O") == 0) {
      memcpy(r->xo, at->x, sizeof(at->x[0])*3);
      r->flags |= HAS_O;
    } else if (strcmp(at->atnm, "CB") == 0) {
      memcpy(r->xcb, at->x, sizeof(at->x[0])*3);
      r->flags |= HAS_CB;
    }
  }
  /* checking integrity */
  for (i = 0; i < c->nres; i++) {
    r = c->res + i;
    if ((r->flags & HAS_BB) != HAS_BB) {
      if (verbose)
        fprintf(stderr, "%s: coordinates of residue %d/%d (%s) is incomplete (0x%X).\n",
          c->file, i+1, c->nres, pdbaaname(r->aa), r->flags);
      c->res[i].broken = 1;
      if (i < c->nres - 1)
        c->res[i+1].broken = 1;
    }
  }
  /* check bond-length */
  for (i = 0; i < c->nres; i++) {
    real x;
    r = c->res + i;
    if (r->broken) continue;
    if (i > 0) {
      x = rv3_dist(c->res[i-1].xc, r->xn);
      if (x < .3 || x > 2.3) {
        if (verbose) {
          const char *aap = pdbaaname(c->res[i-1].aa), *aa = pdbaaname(r->aa);
          fprintf(stderr, "%s: C-N bond between %d (%s) and %d (%s) is broken %g, insert break\n",
            c->file, i, aap, i+1, aa, x);
        }
        c->res[i].broken = 1;
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


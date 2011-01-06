#ifndef PDB_C__
#define PDB_C__

#include "pdb.h"

pdbmodel_t *pdbload0(const char *fname, int verbose)
{
  const int bsiz = 256;
  FILE *fp;
  pdbmodel_t *m;
  pdbatom_t *at;
  int i, j, ir, iro, nline;
  char s[256], resnm[8];

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
  nline = 0;
  ir = -1;
  while (fgets(s, sizeof s, fp)) {
    nline++;
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
    sscanf(s+30, "%lf%lf%lf", at->x, at->x+1, at->x+2);
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

static const char pdb_aanames_[21][4] = {"---",
  "GLY", "ALA", "VAL", "LEU", "ILE", "PRO",
  "THR", "SER", "CYS", "MET", "PHE", "TYR", "TRP",
  "GLU", "GLN", "ASP", "ASN", "ARG", "LYS", "HIS"};

int pdbaaidx(const char *res)
{
  int i;
  for (i = 1; i <= 20; i++)
    if (strcmp(res, pdb_aanames_[i]) == 0) 
      return i;
  return -1;
}

const char *pdbaaname(int i)
{
  if (i <= 0 || i > 20) {
    fprintf(stderr, "invalid pdb id %d\n", i);
    exit(1);
  }
  return pdb_aanames_[i];
}

/* get amino acid backbone */
pdbaabb_t *pdbgetaabb(pdbmodel_t *m, int verbose)
{
  pdbatom_t *at;
  pdbaabb_t *bb;
  pdbaares_t *r;
  int i;

  if ((bb = calloc(1, sizeof(*bb))) == NULL) {
    fprintf(stderr, "no memory for aabb\n");
    return NULL;
  }
  bb->nres = m->nres;
  if ((bb->res = calloc(bb->nres, sizeof(*(bb->res)))) == NULL) {
    fprintf(stderr, "no memory for residues\n");
    free(bb);
    return NULL;
  }
  bb->file = m->file;
  for (i = 0; i < m->n; i++) {
    at = m->at + i;
    r = bb->res + at->rid;
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
  for (i = 0; i < bb->nres; i++) {
    r = bb->res + i;
    if ((r->flags & HAS_BB) != HAS_BB) {
      if (verbose)
        fprintf(stderr, "%s: coordinates of residue %d/%d (%s) is incomplete (0x%X).\n",
          bb->file, i+1, bb->nres, pdbaaname(r->aa), r->flags);
      bb->res[i].broken = 1;
      if (i < bb->nres - 1)
        bb->res[i+1].broken = 1;
    }
  }
  /* check bond-length */
  for (i = 0; i < bb->nres; i++) {
    double x;
#define sqr_(dx) ((dx)*(dx))
#define dist_(a, b) sqrt(sqr_(a[0]-b[0]) + sqr_(a[1]-b[1]) + sqr_(a[2]-b[2]))
    r = bb->res + i;
    if (r->broken) continue;
    if (i > 0) {
      x = dist_(bb->res[i-1].xc, r->xn);
      if (x < .3 || x > 2.3) {
        if (verbose)
          fprintf(stderr, "%s: C-N bond between %d (%s) and %d (%s) is broken %g, insert break\n",
            bb->file, i, pdbaaname(bb->res[i-1].aa), i+1, pdbaaname(r->aa), x);
        bb->res[i].broken = 1;
      }
    }
    x = dist_(r->xn, r->xca);
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: N-CA bond of residue %d (%s) is broken %g\n",
          bb->file, i+1, pdbaaname(r->aa), x);
        fprintf(stderr, "N : %8.3f %8.3f %8.3f\n", r->xn[0], r->xn[1], r->xn[2]);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
      }
      goto ERR;
    }
    x = dist_(r->xca, r->xc);
    if (x < .4 || x > 3.0) {
      if (verbose) {
        fprintf(stderr, "%s: CA-C bond of residue %d (%s) is broken %g\n",
          bb->file, i+1, pdbaaname(r->aa), x);
        fprintf(stderr, "CA: %8.3f %8.3f %8.3f\n", r->xca[0], r->xca[1], r->xca[2]);
        fprintf(stderr, "C : %8.3f %8.3f %8.3f\n", r->xc[0], r->xc[1], r->xc[2]);
      }
      goto ERR;
    }
#undef dist_
#undef sqr_
  }
  if (verbose >= 3) {
    for (i = 0; i < bb->nres; i++)
      printf("%4d: %s\n", i+1, pdbaaname(bb->res[i].aa));
  }
  return bb;
ERR:
  free(bb->res);
  free(bb);
  return NULL;
}

#endif


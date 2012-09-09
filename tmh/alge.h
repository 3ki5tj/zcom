typedef struct {
  int emin, emax, edel, n, order;
  double *dh; /* dH/dE */
  double *cnt; /* histogram */
  double a0, ac; /* alpha */
} algei_t;

static algei_t *algei_open(int emin, int emax, int edel, double a0, double ac,
    double bet0, int order)
{
  algei_t *al;
  int i;
  
  xnew(al, 1);
  al->emin = emin;
  al->emax = emax;
  al->edel = edel; /* cell size */
  al->n = (emax - emin)/edel;
  al->order = order;
  xnew(al->dh, al->n + 1);
  xnew(al->cnt, al->n + 1);
  for (i = 0; i <= al->n; i++) {
    al->dh[i] = bet0;
    al->cnt[i] = 0.;
  }
  al->a0 = a0;
  al->ac = ac;
  return al;
}

static void algei_close(algei_t *al)
{
  free(al->dh);
  free(al->cnt);
  free(al);
}

/* save dh/de file */
static int algei_save(algei_t *al, const char *fn)
{
  FILE *fp;
  int i;
  double lng = 0.;
  
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d %d\n", al->n, al->emin, al->edel);
  for (i = 0; i < al->n; i++) {
    if (i > 0) {
      if (al->order) lng += al->dh[i-1] * al->edel;
      else lng += (al->dh[i-1] + al->dh[i]) * .5 * al->edel;
    }
    fprintf(fp, "%d %.6f %.6f\n", al->emin + i * al->edel, al->dh[i], lng);
  }
  fclose(fp);
  return 0;
}

/* return the H difference zeroth order */
static double algei_dh0(const algei_t *al, int en, int eo)
{
  int ien, ieo, ie, sgn = 1;
  double dh = 0.;
  
  if (eo == en) return 0.;
  if (en < al->emin) {
    dh += (en - al->emin)*al->dh[0];
    en = al->emin;
  } else if (en >= al->emax) {
    dh += (en - al->emax + 1)*al->dh[al->n - 1];
    en = al->emax - 1;
  }
  if (eo == en) return dh;
  if (eo > en) { /* ensure eo < en */
    sgn = eo, eo = en, en = sgn;
    sgn = -1;
  }
  ieo = (eo - al->emin)/al->edel;
  ien = (en - al->emin)/al->edel;
  if (ieo == ien) { /* within a bin */
    dh += sgn*(en - eo)*al->dh[ieo];
  } else { /* cross several bins */
    dh += sgn*(al->emin + (ieo+1)*al->edel - eo)*al->dh[ieo];
    dh += sgn*(en - al->emin - ien*al->edel)*al->dh[ien];
    for (ie = ieo + 1; ie < ien; ie++) dh += sgn*al->edel*al->dh[ie];
  }
  return dh;
}

/* return dh, assuming en and eo are in the same bin */
INLINE double algei_dh1bin(const algei_t *al, int en, int eo, int i)
{
  double slope;
  if (eo == en) return 0; /* in case en == eo == al->emax */
  slope = (al->dh[i+1] - al->dh[i])/al->edel;
  return (al->dh[i] + slope*((eo + en)*.5 - (al->emin + i*al->edel)))*(en - eo);
}

/* return the H difference first order */
static double algei_dh1(const algei_t *al, int en, int eo)
{
  int ien, ieo, ie, sgn = 1;
  double dh = 0.;
  
  if (eo == en) return 0.;
  if (en < al->emin) {
    dh += (en - al->emin)*al->dh[0];
    en = al->emin;
  } else if (en > al->emax) {
    dh += (en - al->emax)*al->dh[al->n - 1];
    en = al->emax;
  }
  if (eo == en) return dh;
  if (eo > en) { /* ensure eo < en */
    sgn = eo, eo = en, en = sgn;
    sgn = -1;
  }
  ieo = (eo - al->emin)/al->edel;
  ien = (en - al->emin)/al->edel;
  if (ieo == ien) { /* within a bin */
    dh += sgn*algei_dh1bin(al, en, eo, ieo);
  } else { /* cross several bins */
    dh += sgn*algei_dh1bin(al, al->emin + (ieo+1)*al->edel, eo, ieo);
    dh += sgn*algei_dh1bin(al, en, al->emin + ien*al->edel, ien);
    for (ie = ieo + 1; ie < ien; ie++)
      dh += sgn*al->edel*(al->dh[ie] + al->dh[ie+1])*.5;
  }
  return dh;
}



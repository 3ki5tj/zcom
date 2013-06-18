#ifndef LJRDF_H__
#define LJRDF_H__
/* Lennard-Jones system: compute the radial distribution function (RDF)
 * using a specially normalized histogram, cf. hist/testrdf.c */



typedef struct {
  hist_t *rdf; /* histogram for radial distribution function */
  int nfr; /* number of frames in rdf */
  lj_t *lj;
} ljrdf_t;



/* open an ljrdf structure, `rmax' can be 0 */
INLINE ljrdf_t *ljrdf_open(lj_t *lj, double dr, double rmax)
{
  ljrdf_t *ljr;

  xnew(ljr, 1);
  ljr->nfr = 0;
  ljr->lj = lj;
  if (rmax <= 0) rmax = lj->l * .5;
  ljr->rdf = hs_open(1, 0, rmax, dr);
  return ljr;
}



INLINE void ljrdf_close(ljrdf_t *ljr)
{
  hs_close(ljr->rdf);
  free(ljr);
}



/* add pairs to the RDF data */
INLINE int ljrdf_add(ljrdf_t *ljr, unsigned flags)
{
  lj_t *lj = ljr->lj;
  int i, j;
  real rc2, dr2, dx[3];
  double dr;

  rc2 = ljr->lj->l/2;
  rc2 = rc2 * rc2;
  for (i = 0; i < lj->n; i++) {
    for (j = i + 1; j < lj->n; j++) {
      if (lj->d == 2)
        dr2 = lj_pbcdist2_2d(dx, lj->x + 2*i, lj->x + 2*j, lj->l);
      else
        dr2 = lj_pbcdist2_3d(dx, lj->x + 3*i, lj->x + 3*j, lj->l);
      if (dr2 >= rc2) continue;
      dr = sqrt(dr2);
      hs_add(ljr->rdf, &dr, 1.0, flags);
    }
  }
  return ++ljr->nfr; /* number of frames */
}



/* header information in writing rdf */
INLINE int ljrdf_fwheader(FILE *fp, void *pdata)
{
  ljrdf_t *ljr = (ljrdf_t *) pdata;
  fprintf(fp, "RDF %d %d %d %.10e | ",
      ljr->nfr, ljr->lj->d, ljr->lj->n, ljr->lj->l);
  return 0;
}



/* header information in reading rdf */
INLINE int ljrdf_frheader(const char *s, void *pdata)
{
  ljrdf_t *ljr = (ljrdf_t *) pdata;
  lj_t *lj = ljr->lj;
  int ret, d, n;
  double l;

  ret = sscanf(s, " RDF %d%d%d%lf | ", &(ljr->nfr), &d, &n, &l);
  die_if (d != lj->d, "dimension mismatch %d vs. %d (file)\n", lj->d, d);
  die_if (n != lj->n, "# of particle mismatch %d vs. %d (file)\n", lj->n, n);
  die_if (fabs(l - lj->l) > 1e-3, "box size mismatch %d vs. %d (file)\n", lj->l, l);
  return (ret == 4) ? 0 : 1;
}



/* normalization */
INLINE double ljrdf_norm(int row, int i, double xmin, double dx, void *pdata)
{
  int npr;
  double x, vsph;
  ljrdf_t *ljr = (ljrdf_t *) pdata;
  lj_t *lj = ljr->lj;

  (void) row;
  x = xmin + i * dx;
  if (lj->d == 2) vsph = 2. * M_PI * dx * (2*x + dx);
  else vsph = 4. * M_PI * dx * (x*(x + dx) + dx*dx/3.);
  npr = lj->n * (lj->n - 1)/2;
  return lj->vol / (vsph * npr * ljr->nfr);
}



/* save rdf, flags can have HIST_NOZEROES */
INLINE int ljrdf_save(ljrdf_t *ljr, const char *fn, unsigned flags)
{
  return hs_savex(ljr->rdf, fn, ljrdf_fwheader, ljrdf_norm, ljr, flags);
}



/* load rdf, flags can have HIST_ADDITION and/or HIST_VERBOSE */
INLINE int ljrdf_load(ljrdf_t *ljr, const char *fn, unsigned flags)
{
  return hs_loadx(ljr->rdf, fn, ljrdf_frheader, ljrdf_norm, ljr, flags);
}



#endif



#ifndef HIST_C__
#define HIST_C__

#include "hist.h"

#define nalloc_(arr, n) \
  if ((arr = calloc((n), sizeof(*arr))) == NULL) { \
    fprintf(stderr, "no memory for %s n = %d\n", #arr, (n)); \
    exit(1); }

static double *gethistsums_(const double *h, int rows, int xcnt)
{
  double *sums;
  int i, j;
  
  nalloc_(sums, rows);
  for (j = 0; j < rows; j++) 
    for (sums[j] = 0., i = 0; i < xcnt; i++)
      sums[j] += h[j*xcnt + i];
  return sums;
}

/* write histograms to file
 * histogram 'h' contains 'rows' histograms,
 * each contains 'xcnt' entries, from 'xmin' to 'xmin+dx*xcnt'
 * (*fwheader) is function to print additional information
 * (*fnorm) is advanced normalization function */
int histsavex(const double *h, int rows, int xcnt, double xmin, double dx,
    unsigned flags, 
    int (*fwheader)(FILE *fp, void *data), 
    double (*fnorm)(int j, int ix, double xmin, double dx, void *data),
    void *pdata,
    const char *fname)
{
  const int version = 0;
  const char *filename;
  FILE *fp;
  int i, j, imax, imin;
  const double *p;
  double *sums, fac, delta;

  filename = (fname != NULL) ? fname : "HIST";

  if ((fp = fopen(filename, "w")) == NULL) {
    printf("cannot write history file [%s].\n", filename);
    return 1;
  }
  
  sums = gethistsums_(h, rows, xcnt);
  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g | ", 
      version, flags, rows, xcnt, xmin, dx);
  for (j = 0; j < rows; j++)
    fprintf(fp, "%g ", sums[j]);
  fprintf(fp, "| ");
  if (fwheader != NULL) (*fwheader)(fp, pdata);
  fprintf(fp, "\n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (j = 0; j < rows; j++) {
    p = h+j*xcnt;

    if (flags & HIST_KEEPRIGHT) {
      imax = xcnt;
    } else { /* trim the right edge */
      for (i = xcnt-1; i >= 0; i--)
        if (p[i] > 0)
          break;
      imax = i+1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge */
      for (i = 0; i < imax; i++)
        if (p[i] > 0)
          break;
      imin = i;
    }

    if (fabs(sums[j]) < 1e-6)
      fac = 1.;
    else fac = 1.0/(sums[j]*dx);

    for (i = imin; i < imax; i++) {
      if ((flags & HIST_NOZEROES) && p[i] < 1e-6)
        continue;
      fprintf(fp,"%g ", xmin+(i+delta)*dx);
      if (flags & HIST_KEEPHIST) 
        fprintf(fp, "%20.14E ", p[i]);
      if (fnorm != NULL) /* advanced normalization */
        fac = (*fnorm)(j, i, xmin, dx, pdata);
      fprintf(fp,"%20.14E %d\n", p[i]*fac, j);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  free(sums);
  return 0;
}

/* load a previous histogram
 * (*frheader) function to read additional header info.
 * (*fnorm) normalization factor */
int histloadx(double *hist, int rows, int xcnt, double xmin, double dx,
    unsigned flags, 
    int (*frheader)(const char *s, void *data), 
    double (*fnorm)(int j, int ix, double xmin, double dx, void *data),
    void *pdata,
    const char *fn)
{
  FILE *fp;
  static char s[4096], *p;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, i1, j, j1, nlin = 0;
  unsigned fflags;
  double x, y, y2, fac, delta, *arr, *sums = NULL;
    
  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;

  if (6 != sscanf(s, " # %d 0x %X | %d%d%lf%lf | %n", &ver, &fflags, &j, &i, &y, &x, &next)
      || i < xcnt || j != rows || fabs(x - dx) > 1e-5) {
    fprintf(stderr, "Error: bins = %d, %d, ng = %d, %d; dx = %g, %g\n", 
        i, xcnt, j, rows, x, dx);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST_KEEPHIST);
  /* scan sums */
  nalloc_(sums, rows);
  for (p = s+next, j = 0; j < rows; j++) {
    if (1 != sscanf(p, "%lf%n", sums + j, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", j, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  /* read the final | */
  next = -1;
  sscanf(p, " | %n", &next);
  if (next == -1) {
    fprintf(stderr, "line %d: cannot read finishing |, s:%s\np:%s\n", j, s, p);
    goto EXIT;
  }
  p += next;
  if (frheader != NULL) {
    if (0 != frheader(p, pdata))
      goto EXIT;
  }


  if (!add) { /* clear histogram */
    for (i = 0; i < rows*xcnt; i++) hist[i] = 0.;
  }
  
  /* loop over j = 0..rows-1 */
  for (j = 0; j < rows; j++) {
    arr = hist + j*xcnt;
    fac = sums[j]*dx;
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; isspace((unsigned char)(*p)) && p >= s; p--) 
        *p = '\0'; /* trim ending */
      if (s[0] == '#' || s[0] == '\0') break;
      if (hashist) {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &y2, &j1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else { /* simple */
        if (3 != sscanf(s, "%lf%lf%d", &x, &y2, &j1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if (j1 != j) {
        fprintf(stderr, "wrong column index %d vs. %d on line %d, s=[%s]\n",
            j1, j, nlin, s);
        goto EXIT;
      }
      i1 = (int)((x - xmin)/dx - delta + .5);
      if (i1 < 0 || i1 >= xcnt) {
        fprintf(stderr, "cannot find index for x = %g, delta = %g, i = %d/%d, on line %d\n", 
            x, delta, i1, xcnt, nlin);
        goto EXIT;
      }
      if (!hashist) {
        if (fnorm != NULL) {
          fac = (*fnorm)(j, i1, xmin, dx, pdata);
          fac = ((fabs(fac) < 1e-8) ? 1. : 1./fac);
        }
        y = y2*fac;
      }
      if (add) arr[i1] += y;
      else arr[i1] = y;
    }
  }
  if (verbose)
    fprintf(stderr, "histogram loaded successfully from %s\n", fn);
  
  if (sums) free(sums);
  return 0;
EXIT:
  if (sums) free(sums);
  /* we always clear histogram on error */
  for (i = 0; i < rows*xcnt; i++) hist[i] = 0.;
  return -1;
}

/* add x of weight w, into histogram h
 * return number of success */
int histadd(const double *xarr, double w, double *h, int rows, 
    int xcnt, double xmin, double dx, unsigned flags)
{
  int j, ix, good = 0, verbose = flags & HIST_VERBOSE;
  double x;

  for (j = 0; j < rows; j++) {
    x = xarr[j];
    if (x < xmin) {
      if (verbose)
       fprintf(stderr, "histadd underflows %d: %g < %g\n", j, x, xmin);
      continue;
    }
    ix = (int)((x - xmin)/dx);
    if (ix >= xcnt) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g > %g\n", j, x, xmin+dx*xcnt);
      continue;
    }
    h[j*xcnt + ix] += w;
    good++;
  }
  return good;
}

/* OO wrappers */
hist_t *hs_initx(int rows, double xmin, double xmax, double dx,
    int (*fwh)(FILE *, void *), int (*frh)(const char*, void *),
    double (*fnorm)(int, int, double, double, void *))
{
  hist_t *hs;

  nalloc_(hs, 1);
  hs->rows = rows;
  hs->xmin = xmin;
  hs->dx   = dx;
  hs->xcnt = (int)((xmax - xmin)/dx + 0.99999999);
  nalloc_(hs->arr, hs->xcnt*hs->rows);
  hs->fwheader = fwh;
  hs->frheader = frh;
  hs->fnorm = fnorm;
  return hs;
}  

void hs_free(hist_t *hs)
{
  free(hs->arr);
  memset(hs, 0, sizeof(*hs));
  free(hs);
}

int hs_savex(const hist_t *hs, const char *fname, void *pdata, unsigned flags)
{
  return histsavex(hs->arr, hs->rows, hs->xcnt, hs->xmin, hs->dx, flags, 
      hs->fwheader, hs->fnorm, pdata, fname);
}

int hs_loadx(hist_t *hs, const char *fname, void *pdata, unsigned flags)
{
  return histloadx(hs->arr, hs->rows, hs->xcnt, hs->xmin, hs->dx, flags,
      hs->frheader, hs->fnorm, pdata, fname);
}

int hs_add(hist_t *hs, const double *x, double w, unsigned flags)
{
  return histadd(x, w, hs->arr, hs->rows, hs->xcnt, hs->xmin, hs->dx, flags);
}

int hs_add1(hist_t *hs, int j, double x, double w, unsigned flags)
{
  return histadd(&x, w, hs->arr + j*hs->xcnt, 1, hs->xcnt, hs->xmin, hs->dx, flags);
}


#undef nalloc_
#endif /* ZCOM_HIST__ */


#ifndef HIST_C__
#define HIST_C__

#include "hist.h"

#ifndef xnew
#define xnew(x, n) \
  if ((n) <= 0) { \
    fprintf(stderr, "cannot allocate %d objects for %s\n", (int) (n), #x); \
    exit(1); \
  } else if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); } 
#endif

/* compute sum, average and standard deviation*/
static double *gethistsums_(const double *h, int rows, int n, 
    double xmin, double dx)
{
  double *sums, *xav, *xdv, x, w;
  int i, r;
  
  xnew(sums, 3*rows);
  xav = sums + rows;
  xdv = xav  + rows;
  for (r = 0; r < rows; r++) {
    sums[r] = xav[r] = xdv[r] = 0.; 
    for (i = 0; i < n; i++) {
      x = xmin + (i+.5)*dx;
      w = h[r*n + i];
      sums[r] += w;
      xav[r]  += w*x;
      xdv[r]  += w*x*x;
    }
    if (sums[r] > 1e-5) {
      xav[r] /= sums[r];
      xdv[r] = sqrt(xdv[r]/sums[r] - xav[r]*xav[r]);
    }
  }
  return sums;
}

static double *gethist2sums_(const double *h, int rows, int n, 
    double xmin, double dx)
{
  double *sums, *xav, *yav, *xdv, *ydv, x, y, w;
  int i, j, r;
  
  xnew(sums, 5*rows);
  xav = sums + rows;
  xdv = sums + rows*2;
  yav = sums + rows*3;
  ydv = sums + rows*4;
  for (r = 0; r < rows; r++) {
    sums[r] = xav[r] = xdv[r] = yav[r] = ydv[r] = 0.;
    for (i = 0; i < n; i++) {
      x = xmin + (i+.5)*dx;
      for (j = 0; j < n; j++) {
        y = xmin + (j+.5)*dx;
        w = h[r*n*n + i*n + j];
        sums[r] += w;
        xav[r]  += w*x;
        xdv[r]  += w*x*x;
        yav[r]  += w*y;
        ydv[r]  += w*y*y;
      }
    }
    if (sums[r] > 1e-5) {
      xav[r] /= sums[r];
      xdv[r] = sqrt(xdv[r]/sums[r] - xav[r]*xav[r]);
      yav[r] /= sums[r];
      ydv[r] = sqrt(ydv[r]/sums[r] - yav[r]*yav[r]);
    }
  }
  return sums;
}


/* write histograms to file
 * histogram 'h' contains 'rows' histograms,
 * each contains 'n' entries, from 'xmin' to 'xmin+dx*n'
 * (*fwheader) is function to print additional information
 * (*fnorm) is advanced normalization function */
int histsavex(const double *h, int rows, int n, double xmin, double dx,
    unsigned flags, 
    int (*fwheader)(FILE *fp, void *data), 
    double (*fnorm)(int r, int ix, double xmin, double dx, void *data),
    void *pdata,
    const char *fn)
{
  const int version = 0;
  const char *filename;
  FILE *fp;
  int i, r, rp, rowp, imax, imin;
  const double *p;
  double sm, *sums, fac, delta;
  double *smtot, *htot = NULL;

  filename = (fn != NULL) ? fn : "HIST";

  if ((fp = fopen(filename, "w")) == NULL) {
    printf("cannot write history file [%s].\n", filename);
    return 1;
  }
  
  sums = gethistsums_(h, rows, n, xmin, dx);
 
  /* compute the overall histogram */
  if (flags & HIST_OVERALL) {
    xnew(htot, n);
    for (i = 0; i < n; i++) htot[i] = 0.;
    
    for (r = 0; r < rows; r++) 
      for (i = 0; i < n; i++)
        htot[i] += h[r*n + i];
    smtot = gethistsums_(htot, 1, n, xmin, dx);
    rowp = rows + 1;
  } else {
    rowp = rows;
  }

  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g | ", 
      version, flags, rows, n, xmin, dx);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r]);
  fprintf(fp, "| ");
  for (r = 0; r < rows; r++) /* average, standard deviation */
    fprintf(fp, "%g %g ", sums[r+rows], sums[r+rows*2]);
  fprintf(fp, "| ");
  if (fwheader != NULL) (*fwheader)(fp, pdata);
  fprintf(fp, "\n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rowp; r++) {
    p = (r == rows) ? htot : (h+r*n);

    if (flags & HIST_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge */
      for (i = n-1; i >= 0; i--)
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

    sm = (r == rows) ? smtot[0] : sums[r];
    if (fabs(sm) < 1e-6) fac = 1.;
    else fac = 1.0/(sm*dx);

    for (i = imin; i < imax; i++) {
      if ((flags & HIST_NOZEROES) && p[i] < 1e-6)
        continue;
      fprintf(fp,"%g ", xmin+(i+delta)*dx);
      if (flags & HIST_KEEPHIST) 
        fprintf(fp, "%20.14E ", p[i]);
      rp = (r == rows) ? (-1) : r;
      if (fnorm != NULL) /* advanced normalization, note the r = -1 case */
        fac = (*fnorm)(rp, i, xmin, dx, pdata);
      fprintf(fp,"%20.14E %d\n", p[i]*fac, rp);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  if (flags & HIST_VERBOSE) {
    fprintf(stderr, "successful wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f av: %10.4f(%10.4f)\n", 
          r, sums[r], sums[r+rows], sums[r+rows*2]);
  }
  free(sums);
  if (flags & HIST_OVERALL) {
    free(htot);
    free(smtot);
  }
  return 0;
}

/* skip a | */
static char *skipabar_(char *p)
{
  int next = -1;
  sscanf(p, " | %n", &next);
  return (next < 0) ? NULL : (p + next);
}

/* load a previous histogram
 * (*frheader) function to read additional header info.
 * (*fnorm) normalization factor */
int histloadx(double *hist, int rows, int n, double xmin, double dx,
    unsigned flags, 
    int (*frheader)(const char *s, void *data), 
    double (*fnorm)(int r, int ix, double xmin, double dx, void *data),
    void *pdata,
    const char *fn)
{
  FILE *fp;
  static char s[40960] = "", *p;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, i1, r, r1, nlin = 0;
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
  if (6 != sscanf(s, " # %d 0x %X | %d%d%lf%lf | %n", &ver, &fflags, &r, &i, &y, &x, &next)
      || i < n || r != rows || fabs(x - dx) > 1e-5) {
    fprintf(stderr, "Error: bins = %d, %d, ng = %d, %d; dx = %g, %g\n", 
        i, n, r, rows, x, dx);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST_KEEPHIST);
  /* scan sums */
  xnew(sums, rows);
  for (p = s+next, r = 0; r < rows; r++) {
    if (1 != sscanf(p, "%lf%n", sums + r, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;
  for (r = 0; r < rows; r++) {
    if (2 != sscanf(p, "%lf%lf%n", &y, &y2, &next)) {
      fprintf(stderr, "cannot read average/stddev from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;
  if (frheader != NULL) {
    if (0 != frheader(p, pdata))
      goto EXIT;
  }


  if (!add) { /* clear histogram */
    for (i = 0; i < rows*n; i++) hist[i] = 0.;
  }
  
  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    arr = hist + r*n;
    fac = sums[r]*dx;
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; isspace((unsigned char)(*p)) && p >= s; p--) 
        *p = '\0'; /* trim ending */
      if (s[0] == '#' || s[0] == '\0') break;
      if (hashist) {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &y2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else { /* simple */
        if (3 != sscanf(s, "%lf%lf%d", &x, &y2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if (r1 < 0) break; /* overall histogram */

      if (r1 < r) {
        fprintf(stderr, "wrong column index %d vs. %d on line %d, s=[%s]\n",
            r1, r, nlin, s);
        goto EXIT;
      } else if (r1 > r) {
        //fprintf(stderr, "jump from column index %d --> %d\n", r, r1);
        r = r1;
        arr = hist + r*n;
        fac = sums[r]*dx;
      }
      i1 = (int)((x - xmin)/dx - delta + .5);
      if (i1 < 0 || i1 >= n) {
        fprintf(stderr, "cannot find index for x = %g, delta = %g, i = %d/%d, on line %d\n", 
            x, delta, i1, n, nlin);
        goto EXIT;
      }
      if (!hashist) {
        if (fnorm != NULL) {
          fac = (*fnorm)(r, i1, xmin, dx, pdata);
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
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  if (sums) free(sums);
  /* we always clear histogram on error */
  for (i = 0; i < rows*n; i++) hist[i] = 0.;
  return -1;
}

/* add x of weight w, into histogram h
 * return number of success */
int histadd(const double *xarr, double w, double *h, int rows, 
    int n, double xmin, double dx, unsigned flags)
{
  int r, ix, good = 0, verbose = flags & HIST_VERBOSE;
  double x;

  for (r = 0; r < rows; r++) {
    x = xarr[r];
    if (x < xmin) {
      if (verbose)
       fprintf(stderr, "histadd underflows %d: %g < %g\n", r, x, xmin);
      continue;
    }
    ix = (int)((x - xmin)/dx);
    if (ix >= n) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g > %g\n", r, x, xmin+dx*n);
      continue;
    }
    h[r*n + ix] += w;
    good++;
  }
  return good;
}

/* OO wrappers */
static void hs_check(const hist_t *hs)
{
  if (hs == NULL) {
    fprintf(stderr, "hist is NULL\n");
    exit(1);
  }
  if (hs->arr == NULL || hs->rows == 0 || hs->n == 0) {
    fprintf(stderr, "hist: arr %p rows %d n %d\n", (void *)(hs->arr), hs->rows, hs->n);
    exit(1);
  }
}

hist_t *hs_initx(int rows, double xmin, double xmax, double dx,
    int (*fwh)(FILE *, void *), int (*frh)(const char*, void *),
    double (*fnorm)(int, int, double, double, void *))
{
  hist_t *hs;

  xnew(hs, 1);
  hs->rows = rows;
  hs->xmin = xmin;
  hs->dx   = dx;
  hs->n = (int)((xmax - xmin)/dx + 0.99999999);
  xnew(hs->arr, hs->n*hs->rows);
  hs->fwheader = fwh;
  hs->frheader = frh;
  hs->fnorm = fnorm;
  return hs;
}  

void hs_free(hist_t *hs)
{
  if (hs) {
    if (hs->arr) free(hs->arr);
    memset(hs, 0, sizeof(*hs));
    free(hs);
  }
}

int hs_savex(const hist_t *hs, const char *fn, void *pdata, unsigned flags)
{
  hs_check(hs);
  return histsavex(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags, 
      hs->fwheader, hs->fnorm, pdata, fn);
}

int hs_loadx(hist_t *hs, const char *fn, void *pdata, unsigned flags)
{
  hs_check(hs);
  return histloadx(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags,
      hs->frheader, hs->fnorm, pdata, fn);
}

int hs_add(hist_t *hs, const double *x, double w, unsigned flags)
{
  hs_check(hs);
  return histadd(x, w, hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags);
}

int hs_add1(hist_t *hs, int r, double x, double w, unsigned flags)
{
  hs_check(hs);
  if (r >= hs->rows || r < 0) {
    fprintf(stderr, "bad row index %d\n", r);
    exit(1);
  }
  return histadd(&x, w, hs->arr + r*hs->n, 1, hs->n, hs->xmin, hs->dx, flags);
}


/* write 'rows' 2d n^2 histograms to file */
int hist2save(const double *h, int rows, int n, double xmin, double dx,
    unsigned flags, const char *fn)
{
  const int version = 0;
  const char *filename;
  FILE *fp;
  int i, j, r, imax, imin, jmax, jmin, n2;
  const double *p;
  double *sums, fac, delta;

  filename = (fn != NULL) ? fn : "HIST2";

  if ((fp = fopen(filename, "w")) == NULL) {
    printf("cannot write history file [%s].\n", filename);
    return 1;
  }
  
  n2 = n*n;
  sums = gethist2sums_(h, rows, n, xmin, dx);
  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g | ", 
      version, flags, rows, n, xmin, dx);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r]);
  fprintf(fp, " | ");
  for (r = 0; r < rows; r++) /* averages and standard deviations */
    fprintf(fp, "%g %g %g %g ", sums[r+rows], sums[r+rows*2],
        sums[r+rows*3], sums[r+rows*4]);
  fprintf(fp, "| \n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rows; r++) {
    p = h+r*n2;

    if (flags & HIST_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge of i */
      for (i = n-1; i >= 0; i--) {
        for (j = 0; j < n; j++)
          if (p[i*n + j] > 0) break;
        if (j < n) break;
      }
      imax = i+1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge of i */
      for (i = 0; i < imax; i++) {
        for (j = 0; j < n; j++)
          if (p[i*n + j] > 0) break;
        if (j < n) break;
      }
      imin = i;
    }

    if (flags & HIST_KEEPRIGHT2) {
      jmax = n;
    } else { /* trim the right edge of j */
      for (j = n-1; j >= 0; j--) {
        for (i = imin; i < imax; i++)
          if (p[i*n + j] > 0) break;
        if (i < imax) break;
      }
      jmax = j+1;
    }
    
    if (flags & HIST_KEEPLEFT2) {
      jmin = 0;
    } else { /* trim the left edge of j */
      for (j = 0; j < jmax; j++) {
        for (i = imin; i < imax; i++)
          if (p[i*n + j] > 0) break;
        if (i < imax) break;
      }
      jmin = j;
    }

    if (fabs(sums[r]) < 1e-6) fac = 1.;
    else fac = 1.0/(sums[r]*dx*dx);

    for (i = imin; i < imax; i++) { 
      for (j = jmin; j < jmax; j++) {
        double x, y;
        if ((flags & HIST_NOZEROES) && p[i] < 1e-6)
          continue;
        x = xmin + (i+delta)*dx;
        y = xmin + (j+delta)*dx;
        fprintf(fp,"%g %g ", x, y);
        if (flags & HIST_KEEPHIST) 
          fprintf(fp, "%20.14E ", p[i*n+j]);
        fprintf(fp,"%20.14E %d\n", p[i*n+j]*fac, r);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp, "\n#\n");
  }
  fclose(fp);
  if (flags & HIST_VERBOSE) {
    fprintf(stderr, "successful wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f xav: %10.4f(%10.4f) yav: %10.4f(%10.4f)\n", 
          r, sums[r], sums[r+rows], sums[r+rows*2], sums[r+rows*3], sums[r+rows*4]);
  }
  free(sums);
  return 0;
}

int hist2load(double *hist, int rows, int n, double xmin, double dx,
    unsigned flags, const char *fn)
{
  FILE *fp;
  static char s[40960] = "", *p;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, j, r, r1, n2, nlin = 0;
  unsigned fflags;
  double x, y, g, g2, fac, delta, *arr, *sums = NULL;
    
  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  
  n2 = n*n;
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;
  if (6 != sscanf(s, " # %d 0x %X | %d%d%lf%lf | %n", &ver, &fflags, &r, &i, &y, &x, &next)
      || i < n || r != rows || fabs(x - dx) > 1e-5) {
    fprintf(stderr, "Error: bins = %d, %d, ng = %d, %d; dx = %g, %g\n", 
        i, n, r, rows, x, dx);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST_KEEPHIST);
  /* scan sums */
  xnew(sums, rows);
  for (p = s+next, r = 0; r < rows; r++) {
    if (1 != sscanf(p, "%lf%n", sums + r, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;
  for (r = 0; r < rows; r++) {
    if (4 != sscanf(p, "%lf%lf%lf%lf%n", &x, &y, &g, &g2, &next)) {
      fprintf(stderr, "cannot read average/stddev from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = skipabar_(p)) == NULL) goto EXIT;

  if (!add) { /* clear histogram */
    for (i = 0; i < rows*n2; i++) hist[i] = 0.;
  }
  
  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    arr = hist + r*n2;
    fac = sums[r]*(dx*dx);
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; isspace((unsigned char)(*p)) && p >= s; p--) 
        *p = '\0'; /* trim ending */
      if (s[0] == '#') break;
      if (s[0] == '\0') continue;
      
      if (hashist) {
        if (5 != sscanf(s, "%lf%lf%lf%lf%d", &x, &y, &g, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if (r1 != r) {
        fprintf(stderr, "wrong column index %d vs. %d on line %d\n",
          r1, r, nlin);
        goto EXIT;
      }
      i = (int)((x - xmin)/dx - delta + .5);
      if (i < 0 || i >= n) {
        fprintf(stderr, "cannot find index for x = %g\n", x);
        goto EXIT;
      }
      j = (int)((y - xmin)/dx - delta + .5);
      if (j < 0 || j >= n) {
        fprintf(stderr, "cannot find index for y = %g\n", y);
        return -1;
      }
      if (!hashist) {
        g = g2*fac;
      }
      if (add) arr[i*n+j] += g;
      else arr[i*n+j] = g;
    }
  }
  if (verbose) fprintf(stderr, "%s loaded successfully\n", fn);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);  
  if (sums) free(sums);
  for (i = 0; i < rows*n2; i++) hist[i] = 0.;
  return -1;
}

/* add (xarr[skip*r], yarr[skip*r]) of weight w, into histogram h
 * return number of success */
int hist2add(const double *xarr, const double *yarr, int skip,
    double w, double *h, int rows, 
    int n, double xmin, double dx, unsigned flags)
{
  int r, ix, iy, good = 0, verbose = flags & HIST_VERBOSE;
  double x, y;

  for (r = 0; r < rows; r++) {
    x = xarr[skip*r];
    y = yarr[skip*r];
    if (x < xmin || y < xmin) {
      if (verbose)
       fprintf(stderr, "histadd underflows %d: %g or %g < %g\n", r, x, y, xmin);
      continue;
    }
    ix = (int)((x - xmin)/dx);
    iy = (int)((y - xmin)/dx);
    if (ix >= n || iy >= n) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g or %g > %g\n", r, x, y, xmin+dx*n);
      continue;
    }
    h[r*n*n + ix*n+iy] += w;
    good++;
  }
  return good;
}

static void hs2_check(const hist2_t *hs)
{
  if (hs == NULL) {
    fprintf(stderr, "hist2 is NULL\n");
    exit(1);
  }
  if (hs->arr == NULL || hs->rows == 0 || hs->n == 0) {
    fprintf(stderr, "hist2: arr %p rows %d n %d\n", (void *)(hs->arr), hs->rows, hs->n);
    exit(1);
  }
}

hist2_t *hs2_init(int rows, double xmin, double xmax, double dx)
{
  hist2_t *hs2;

  xnew(hs2, 1);
  hs2->rows = rows;
  hs2->xmin = xmin;
  hs2->dx   = dx;
  hs2->n = (int)((xmax - xmin)/dx + 0.99999999);
  xnew(hs2->arr, hs2->n*hs2->n*hs2->rows);
  return hs2;
}

void hs2_free(hist2_t *hs2)
{
  if (hs2) {
    if (hs2->arr) free(hs2->arr);
    memset(hs2, 0, sizeof(*hs2));
    free(hs2);
  }
}

int hs2_save(const hist2_t *hs, const char *fn, unsigned flags)
{
  hs2_check(hs);
  return hist2save(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, 
      flags, fn);
}

int hs2_load(hist2_t *hs, const char *fn, unsigned flags)
{
  hs2_check(hs);
  return hist2load(hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, 
      flags, fn);
}

int hs2_add(hist2_t *hs, const double *x, const double *y, int skip, double w, unsigned flags)
{
  hs2_check(hs);
  return hist2add(x, y, skip, w, hs->arr, hs->rows, hs->n, hs->xmin, hs->dx, flags);
}

int hs2_add1(hist2_t *hs, int r, double x, double y, double w, unsigned flags)
{
  hs2_check(hs);
  return hist2add(&x, &y, 1, w, hs->arr+r*hs->n*hs->n, 1, hs->n, hs->xmin, hs->dx, flags);
}

#endif /* ZCOM_HIST__ */


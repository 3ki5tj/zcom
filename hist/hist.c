#ifndef HIST_C__
#define HIST_C__
#include "hist.h"

/* compute sum, average and variance */
INLINE double *histgetsums_(const double *h, int rows, int n, 
    double xmin, double dx, double *sums)
{
  double *xav, *xxc, x, w;
  int i, r;
  
  xav = sums + rows;
  xxc = xav  + rows;
  for (r = 0; r < rows; r++) {
    sums[r] = xav[r] = xxc[r] = 0.; 
    for (i = 0; i < n; i++) {
      x = xmin + (i+.5)*dx;
      w = h[r*n + i];
      sums[r] += w;
      xav[r]  += w*x;
      xxc[r]  += w*x*x;
    }
    if (sums[r] > 1e-5) {
      xav[r] /= sums[r];
      xxc[r] = xxc[r]/sums[r] - xav[r]*xav[r];
    }
  }
  return sums;
}


/* get average of a certain `row' */
INLINE double hs_getave(const hist_t *hs, int row, double *sum, double *var)
{
  double arr[3];

  histgetsums_(hs->arr + row * hs->n, 1, hs->n, hs->xmin, hs->dx, arr);
  if (sum) *sum = arr[0];
  if (var) *var = arr[2];
  return arr[1];
}


/* write histograms to file
 * histogram 'h' contains 'rows' histograms,
 * each contains 'n' entries, from 'xmin' to 'xmin+dx*n'
 * (*fwheader) is function to print additional information
 * (*fnorm) is advanced normalization function */
INLINE int histsavex(const double *h, int rows, int n, double xmin, double dx,
    unsigned flags, 
    int (*fwheader)(FILE *fp, void *data), 
    double (*fnorm)(int r, int ix, double xmin, double dx, void *data),
    void *pdata,
    const char *fn)
{
  const int version = 0;
  FILE *fp;
  int i, r, rp, rowp, imax, imin;
  const double *p;
  double sm, *sums, fac, delta;
  double smtot[3], *htot = NULL;

  if (fn == NULL) fn = "HIST";
  xfopen(fp, fn, "w", return -1);
 
  /* get statistics */ 
  xnew(sums, rows * 3);
  histgetsums_(h, rows, n, xmin, dx, sums);
 
  /* compute the overall histogram */
  if (flags & HIST_OVERALL) {
    xnew(htot, n); /* the overall histogram */
    for (i = 0; i < n; i++) htot[i] = 0.;
    
    for (r = 0; r < rows; r++) 
      for (i = 0; i < n; i++)
        htot[i] += h[r*n + i];
    histgetsums_(htot, 1, n, xmin, dx, smtot);
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
    fprintf(fp, "%g %g ", sums[r+rows], sqrt(sums[r+rows*2]));
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
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f av: %10.4f(%10.4f)\n", 
          r, sums[r], sums[r+rows], sums[r+rows*2]);
  }
  free(sums);
  if (flags & HIST_OVERALL) {
    free(htot);
  }
  return 0;
}

/* fetch histogram size */
INLINE int histgetinfo(const char *fn, int *row, double *xmin, double *xmax, double *xdel,
    int *version, unsigned *fflags)
{
  FILE *fp;
  char s[1024];
  int n;
  
  xfopen(fp, fn, "r", return -1);
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);    
    return -1;
  }
  if (6 != sscanf(s, "# %d 0x %X | %d %d %lf %lf ", 
        version, fflags, row, &n, xmin, xdel)) {
    fprintf(stderr, "%s: bad first line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  *xmax = *xmin + *xdel * n;
  fclose(fp);
  return 0;
}

/* skip a | */
INLINE char *skipabar_(char *p)
{
  int next = -1;
  sscanf(p, " | %n", &next);
  return (next < 0) ? NULL : (p + next);
}

/* load a previous histogram
 * (*frheader) function to read additional header info.
 * (*fnorm) normalization factor
 * flags can have HIST_ADDITION and/or HIST_VERBOSE */
INLINE int histloadx(double *hist, int rows, int n, double xmin, double dx,
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
    
  xfopen(fp, fn, "r", return -1);
  
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
  fclose(fp);
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
INLINE int histadd(const double *xarr, double w, double *h, int rows, 
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


INLINE double *hist2getsums_(const double *h, int rows, int n,
    double xmin, double dx, int m, double ymin, double dy, double *sums)
{
  double *xav, *yav, *xxc, *yyc, *xyc, x, y, w;
  int i, j, r;
  
  xav = sums + rows;
  xxc = sums + rows*2;
  yav = sums + rows*3;
  yyc = sums + rows*4;
  xyc = sums + rows*5;
  for (r = 0; r < rows; r++) {
    sums[r] = xav[r] = xxc[r] = yav[r] = yyc[r] = 0.;
    for (i = 0; i < n; i++) {
      x = xmin + (i+.5)*dx;
      for (j = 0; j < m; j++) {
        y = ymin + (j+.5)*dy;
        w = h[r*n*m + i*m + j];
        sums[r] += w;
        xav[r]  += w*x;
        xxc[r]  += w*x*x;
        yav[r]  += w*y;
        yyc[r]  += w*y*y;
        xyc[r]  += w*x*y;
      }
    }
    if (sums[r] > 1e-5) {
      xav[r] /= sums[r];
      xxc[r] = sqrt(xxc[r]/sums[r] - xav[r]*xav[r]);
      yav[r] /= sums[r];
      yyc[r] = sqrt(yyc[r]/sums[r] - yav[r]*yav[r]);
      xyc[r] = xyc[r]/sums[r] - xav[r]*yav[r];
    }
  }
  return sums;
}


/* write 'rows' 2d n^2 histograms to file */
INLINE int hist2save(const double *h, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  const int version = 1; /* v1 allows different dimension in x and y */
  FILE *fp;
  int i, j, r, imax, imin, jmax, jmin, nm;
  const double *p;
  double *sums, fac, delta;

  if (fn == NULL) fn = "HIST2";
  xfopen(fp, fn, "w", return 1);
  
  nm = n*m;
  xnew(sums, rows * 6);
  hist2getsums_(h, rows, n, xmin, dx, m, ymin, dy, sums);
  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g %d %g %g | ", 
      version, flags, rows, n, xmin, dx, m, ymin, dy);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r]);
  fprintf(fp, " | ");
  for (r = 0; r < rows; r++) /* averages and standard deviations */
    fprintf(fp, "%g %g %g %g %g ", sums[r+rows], sums[r+rows*2],
        sums[r+rows*3], sums[r+rows*4], sums[r+rows*5]);
  fprintf(fp, "| \n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rows; r++) { /* the rth data set */
    p = h + r*nm;

    if (flags & HIST_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge of i */
      for (i = n-1; i >= 0; i--) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imax = i+1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge of i */
      for (i = 0; i < imax; i++) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imin = i;
    }

    if (flags & HIST_KEEPRIGHT2) {
      jmax = m;
    } else { /* trim the right edge of j */
      for (j = m-1; j >= 0; j--) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmax = j+1;
    }
    
    if (flags & HIST_KEEPLEFT2) {
      jmin = 0;
    } else { /* trim the left edge of j */
      for (j = 0; j < jmax; j++) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmin = j;
    }

    if (fabs(sums[r]) < 1e-6) fac = 1.;
    else fac = 1.0/(sums[r]*dx*dy);

    for (i = imin; i < imax; i++) { 
      for (j = jmin; j < jmax; j++) {
        double x, y;
        if ((flags & HIST_NOZEROES) && p[i*m + j] < 1e-16)
          continue;
        x = xmin + (i+delta)*dx;
        y = ymin + (j+delta)*dy;
        fprintf(fp,"%g %g ", x, y);
        if (flags & HIST_KEEPHIST) 
          fprintf(fp, "%20.14E ", p[i*m+j]);
        fprintf(fp,"%20.14E %d\n", p[i*m+j]*fac, r);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp, "\n#\n");
  }
  fclose(fp);
  if (flags & HIST_VERBOSE) {
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f xav: %10.4f(%10.4f) yav: %10.4f(%10.4f)\n", 
          r, sums[r], sums[r+rows], sums[r+rows*2], sums[r+rows*3], sums[r+rows*4]);
  }
  free(sums);
  return 0;
}

INLINE int hist2load(double *hist, int rows, int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags, const char *fn)
{
  FILE *fp;
  static char s[40960] = "", *p;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, j, r, r1, nm, nlin = 0;
  unsigned fflags;
  double x, y, g, g2, fac, delta, *arr, *sums = NULL;
  double xmin1, dx1, ymin1, dy1;
    
  xfopen(fp, fn, "r", return -1);
  
  nm = n*m;
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;
  if (9 != sscanf(s, " # %d 0x %X | %d%d%lf%lf%d%lf%lf | %n", &ver, &fflags, &r, 
        &i, &xmin1, &dx1, &j, &ymin1, &dy1, &next)
      || i < n || j < m || r != rows
      || fabs(dx1 - dx) > 1e-5 || fabs(dy1 - dy) > 1e-5 ) {
    fprintf(stderr, "Error: bins %d, %d; %d, %d; ng %d, %d; dx %g, %g; dy %g, %g\n", 
        i, n, j, m, r, rows, dx1, dx, dy1, dy);
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
    for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  }
  
  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    arr = hist + r*nm;
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
      j = (int)((y - ymin)/dy - delta + .5);
      if (j < 0 || j >= m) {
        fprintf(stderr, "cannot find index for y = %g\n", y);
        return -1;
      }
      if (!hashist) {
        g = g2*fac;
      }
      if (add) arr[i*m+j] += g;
      else arr[i*m+j] = g;
    }
  }
  if (verbose) fprintf(stderr, "%s loaded successfully\n", fn);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);  
  if (sums) free(sums);
  for (i = 0; i < rows*nm; i++) hist[i] = 0.;
  return -1;
}

/* add (xarr[skip*r], yarr[skip*r]) of weight w, into histogram h
 * return number of success */
INLINE int hist2add(const double *xarr, const double *yarr, int skip,
    double w, double *h, int rows, 
    int n, double xmin, double dx,
    int m, double ymin, double dy, unsigned flags)
{
  int r, ix, iy, good = 0, verbose = flags & HIST_VERBOSE;
  double x, y;

  for (r = 0; r < rows; r++) {
    x = xarr[skip*r];
    y = yarr[skip*r];
    if (x < xmin || y < ymin) {
      if (verbose)
       fprintf(stderr, "histadd underflows %d: %g < %g or %g < %g\n",
          r, x, xmin, y, ymin);
      continue;
    }
    ix = (int)((x - xmin)/dx);
    iy = (int)((y - ymin)/dy);
    if (ix >= n || iy >= m) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g > %g or %g > %g\n",
            r, x, xmin + dx*n, y, ymin + dy*m);
      continue;
    }
    h[r*n*m + ix*m + iy] += w;
    good++;
  }
  return good;
}


#endif /* ZCOM_HIST__ */

